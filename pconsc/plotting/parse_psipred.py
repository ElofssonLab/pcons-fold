#!/usr/bin/env python


def horizontal(hfile):

    """Reads psipred output .horiz file.
    @param  hfile   psipred .horiz file
    @return secondary structure string.
    """
    
    result = ''
    for line in hfile:
        line_arr = line.strip().split(' ')
        if line_arr[0] == 'Pred:' and len(line_arr) > 1:
            result += line_arr[1]
    return result


def horizontal_conf(hfile):

    """Reads psipred output .horiz file.
    @param  hfile   psipred .horiz file
    @return secondary structure string.
    """
    
    result = ''
    for line in hfile:
        line_arr = line.strip().split(' ')
        if line_arr[0] == 'Conf:' and len(line_arr) > 1:
            result += line_arr[1]
    return result


def horizontal_seq(hfile):

    """Reads psipred output .horiz file.
    @param  hfile   psipred .horiz file
    @return amino acid sequence.
    """

    result = ''
    for line in hfile:
        line_arr = line.strip().split(' ')
        if line_arr[0] == 'AA:':
            result += line_arr[1]
    return result
