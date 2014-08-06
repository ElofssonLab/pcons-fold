def read(rundir, scorefile):

    scores_dict = {}
    
    scorefile.readline()

    for line in scorefile:
        line_arr = line.strip().split()
        tag = line_arr[-1]
        rundir_tag = '%s/%s' % (rundir, tag)
        scores_dict[rundir_tag] = map(float, line_arr[1:-1])

    return scores_dict


def read_successful(rundir, scorefile):

    scores_dict = {}
    
    scorefile.readline()

    for line in scorefile:
        line_arr = line.strip().split()
        tag = line_arr[-1]
        # when reading regular score files
        if tag.startswith('S_'):
            rundir_tag = '%s/%s' % (rundir, tag)
            scores_dict[rundir_tag] = map(float, line_arr[1:-1])
        # in case of rescoring
        elif not tag.startswith('F_'):
            rundir_tag = '%s/%s' % (rundir, tag)
            scores_dict[rundir_tag] = map(float, line_arr[1:-1])

    if len(scores_dict) == 0:
        scorefile.seek(0)
        scores_dict = read(rundir, scorefile)

    return scores_dict
