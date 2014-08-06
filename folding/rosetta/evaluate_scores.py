import sys
import os
import shutil
import operator
from collections import defaultdict
import subprocess

import parse_rosetta_scores
import parse_tmscore
import fix_numbering
from localconfig import *

def check_output(command):
    return subprocess.Popen(command, stdout=subprocess.PIPE).communicate()[0]



def get_best_models(num, rundir_i, scorefile):

    scores_dict = parse_rosetta_scores.read_successful(rundir_i, scorefile)
    scorefile.close()
    scores_sorted = sorted(scores_dict.iteritems(), key=operator.itemgetter(1))

    return scores_sorted[:num]


def get_best_of_all_runs(num, nruns, rundir, scorefile_name='score.fsc'):

    all_scores = []
    
    for i in xrange(nruns):
        rundir_i = '%s/run_%d' % (rundir, i+1)
        scorefile = open('%s/%s' % (rundir_i, scorefile_name), 'r')
        curr_scores = get_best_models(num, rundir_i, scorefile)
        all_scores += curr_scores

    all_scores_sorted = sorted(all_scores, key=lambda(x):x[1][0])

    return all_scores_sorted[:num]


def write_table(scores_sorted, outfile):

    for (rundir_tag, score_list) in scores_sorted:
        score_str = '\t'.join(score_list)
        outfile.write('%s\t%s\n' % (score_str, rundir_tag))

    outfile.close()


def extract_structures(scores_sorted):
    i = 0
    #print scores_sorted
    currdir = os.getcwd()
    for (rundir_tag, score_list) in scores_sorted:
        i += 1
        tag = rundir_tag.split('/')[-1]
        rundir = '/'.join(rundir_tag.split('/')[:-1])
        run = rundir_tag.split('/')[-2]
        os.chdir('%s/' % rundir)
        if os.path.exists('../%d.%s.%s.pdb' % (i, run, tag)):
            os.chdir(currdir)
            continue 
        check_output(['%s/extract_pdbs.linuxgccrelease' % rosetta_binary_dir, 
              '-in:file:silent', 'default.out', 
              '-in:file:tags', '%s' % tag, 
              '-database', rosetta_db_dir])
        shutil.move('%s.pdb' % tag, '../%d.%s.%s.pdb' % (i, run, tag))
        os.chdir(currdir)


def relax_structures(scores_sorted):
    i = 0
    currdir = os.getcwd()
    for (rundir_tag, score_list) in scores_sorted:
        i += 1
        tag = rundir_tag.split('/')[-1]
        rundir = '/'.join(rundir_tag.split('/')[:-1])
        run = rundir_tag.split('/')[-2]
        os.chdir('%s/' % '/'.join(rundir.split('/')[:-1]))
        if os.path.exists('%d.%s.%s_0001.pdb' % (i, run, tag)):
            os.chdir(currdir)
            continue 
        check_output(['%s/relax.linuxgccrelease' % rosetta_binary_dir, 
            '-in:file:s', '%d.%s.%s.pdb' % (i, run, tag), 
            '-in:file:fullatom', 
            '-relax:quick',
            '-database', rosetta_db_dir,
            '-relax:constrain_relax_to_start_coords'])
        os.chdir(currdir)


def rescore_structures(scores_sorted, relax_flag):

    model_lst = []

    i = 1
    for (rundir_tag, score_list) in scores_sorted:
        tag = rundir_tag.split('/')[-1]
        rundir = '/'.join(rundir_tag.split('/')[:-1])
        run = rundir_tag.split('/')[-2]
        os.chdir('%s/%s/' % (rootdir, '/'.join(rundir.split('/')[:-1])))
        if relax_flag:
            model_filename = '%d.%s.%s_0001.pdb' % (i, run, tag)
        else:
            model_filename = '%d.%s.%s.pdb' % (i, run, tag)
        os.chdir('%s/%s/' % (rootdir, '/'.join(rundir.split('/')[:-1])))
        check_output(['%s/score.linuxgccrelease' % rosetta_binary_dir, 
            '-in:file:s', model_filename, 
            '-out:nooutput', 
            '-database', rosetta_db_dir])
        os.chdir(rootdir)
        model_lst.append(model_filename)
        os.chdir('%s' % rootdir)
        i += 1
    

def compare_to_native(scores_sorted, relax_flag, rescore_flag, native_fname='native.pdb', name='', chain='A'):

    scores = defaultdict(list)
    rosetta_scores = []
    if not name:
       name = scores_sorted[0][0].split('/')[0]
    i = 0
    currdir = os.getcwd()
    for (rundir_tag, score_list) in scores_sorted:
        i += 1
        rosetta_scores.append(score_list[0])
        tag = rundir_tag.split('/')[-1]
        if rescore_flag:
            configdir = '/'.join(rundir_tag.split('/')[:-1])
            protdir = '/'.join(rundir_tag.split('/')[:-2])
            run = ''
        else: 
            configdir = '/'.join(rundir_tag.split('/')[:-2])
            protdir = '/'.join(rundir_tag.split('/')[:-3])
            run = rundir_tag.split('/')[-2]
        os.chdir('%s/' % configdir)# (rootdir, '/'.join(rundir.split('/')[:-1])))
        if relax_flag:
            if rescore_flag:
                model_filename = '%s.pdb' % '_'.join(tag.split('_')[:-1])
                score_filename = '%s.rescore.TMscore' % '_'.join(tag.split('_')[:-1])
            else:
                model_filename = '%d.%s.%s_0001.pdb' % (i, run, tag)
                score_filename = '%d.%s.%s_0001.TMscore' % (i, run, tag)
        else:
            model_filename = '%d.%s.%s.pdb' % (i, run, tag)
            score_filename = '%d.%s.%s.TMscore' % (i, run, tag)

        native_aligned_fname = native_fname[:-4] + '.aligned.pdb'
        fix_numbering.fix(model_filename, native_fname, native_aligned_fname, chain2=chain)
        
        score_str = check_output(['%s' % tmscore_binary, model_filename, native_aligned_fname])
        score_f = open(score_filename, 'w')
        score_f.write(score_str)
        score_f.close()

        tmp_scores = parse_tmscore.read(open(score_filename))
        for key, score in tmp_scores.iteritems():
            scores[key].append(score)

        os.chdir(currdir)
    
    avg_scores = {}
    for key, score_lst in scores.iteritems():
        avg_scores[key] = sum(score_lst) / len(score_lst)

    outstr = ''
    outstr += '%s\tRosetta-score\t%s\n' % (name, '\t'.join(map(str, rosetta_scores)))
    outstr += '%s\tRMSD\t%s\n' % (name, '\t'.join(map(str, scores['RMSD'])))
    outstr += '%s\tTM-score\t%s\n' % (name, '\t'.join(map(str, scores['TM-score'])))
    outstr += '%s\tMaxSub\t%s\n' % (name, '\t'.join(map(str, scores['MaxSub'])))
    outstr += '%s\tGDT-TS\t%s\n' % (name, '\t'.join(map(str, scores['GDT-TS'])))
    outstr += '%s\tGDT-HA\t%s\n' % (name, '\t'.join(map(str, scores['GDT-HA'])))
    print outstr
    return outstr


if __name__ == '__main__':

    if len(sys.argv) != 6:
        sys.stderr.write('Usage: evaluate_scores.py <number of models> <number of runs> <run directory> <score table filename> <relax flag>')
    
    num = int(sys.argv[1])
    nruns = int(sys.argv[2])
    #rundir = sys.argv[3]
    top_scores_filename = sys.argv[4]
    relax_flag = bool(int(sys.argv[5]))
    all_scores_sorted = get_best_of_all_runs(num, nruns, rundir)

    #print all_scores_sorted
    top_rundir_tag = all_scores_sorted[0][0]
    top_tag = top_rundir_tag.split('/')[-1]
    top_rundir = '/'.join(top_rundir_tag.split('/')[:2])
    #print top_rundir
    #print top_tag
    #write_table(all_scores_sorted, open(top_scores_filename,'w'))
    extract_structures(all_scores_sorted)
    if relax_flag:
        relax_structures(all_scores_sorted)
    #print 'Name\tRMSD\tTM-score\tMaxSub\tGDT-TS\tGDT-HA'
    compare_to_native(all_scores_sorted, relax_flag, False)
    #rescore_structures(all_scores_sorted, relax_flag)
    #configdir = '/'.join(rundir.split('/')[:-1])
    #all_rescores_sorted = get_best_models(num, configdir, open('%s/default.sc' % configdir, 'r'))
    #compare_to_native(all_rescores_sorted, relax_flag, True)
     



