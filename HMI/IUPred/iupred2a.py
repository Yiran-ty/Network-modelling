#!/usr/bin/python3

import textwrap
import math
import os
from pathlib import Path


def avg(lst):
    return sum(lst) / len(lst)


def aa_freq(_seq):
    _freq = {}
    for _aa in _seq:
        if _aa in _freq:
            _freq[_aa] += 1
        else:
            _freq[_aa] = 1
    for _aa, _ins in _freq.items():
        _freq[_aa] = _ins / len(_seq)
    return _freq


def read_matrix(matrix_file):
    _mtx = {}
    with open(matrix_file, "r") as _fhm:
        for _line in _fhm:
            if _line.split()[0] in _mtx:
                _mtx[_line.split()[0]][_line.split()[1]] = float(_line.split()[2])
            else:
                _mtx[_line.split()[0]] = {}
                _mtx[_line.split()[0]][_line.split()[1]] = float(_line.split()[2])
    return _mtx


def read_histo(histo_file):
    hist = []
    h_min = float("inf")
    h_max = -float("inf")
    with open(histo_file, "r") as fnh:
        for _line in fnh:
            if _line.startswith("#"):
                continue
            if float(_line.split()[1]) < h_min:
                h_min = float(_line.split()[1])
            if float(_line.split()[1]) > h_max:
                h_max = float(_line.split()[1])
            hist.append(float(_line.split()[-1]))
    h_step = (h_max - h_min) / (len(hist))
    return hist, h_min, h_max, h_step


def smooth(energy_list, window):
    weighted_energy_score = [0] * len(energy_list)
    for idx in range(len(energy_list)):
        weighted_energy_score[idx] = avg(energy_list[max(0, idx - window):min(len(energy_list), idx + window + 1)])
    return weighted_energy_score


def read_seq(fasta_file):
    _seq = ""
    if fasta_file.endswith(".fasta"):
        with open(fasta_file) as file_handler:
            for _line in file_handler:
                if _line.startswith(">"):
                    continue
                _seq += _line.strip()

    return _seq


def iupred(folder, seq, mode='short'):
    if mode == "short":
        lc = 1
        uc = 25
        wc = 10
        p = Path(os.path.join(folder, 'iupred2_short_energy_matrix'))
        mtx = read_matrix(p)
        histo, histo_min, histo_max, histo_step = read_histo(os.path.join(folder, 'short_histogram'))

    elif mode == 'glob':
        lc = 1
        uc = 100
        wc = 15
        p = Path(os.path.join(folder, 'iupred2_long_energy_matrix'))
        mtx = read_matrix(p)
        histo, histo_min, histo_max, histo_step = read_histo(os.path.join(folder, 'long_histogram'))

    else:
        lc = 1
        uc = 100
        wc = 10
        p = Path(os.path.join(folder, 'iupred2_long_energy_matrix'))
        p2 = Path(os.path.join(folder, 'long_histogram'))
        mtx = read_matrix(p)
        histo, histo_min, histo_max, histo_step = read_histo(p2)

    unweighted_energy_score = [0] * len(seq)
    weighted_energy_score = [0] * len(seq)
    iupred_score = [0] * len(seq)

    for idx in range(len(seq)):
        freq_dct = aa_freq(seq[max(0, idx - uc):max(0, idx - lc)] + seq[idx + lc + 1:idx + uc + 1])
        for aa, freq in freq_dct.items():
            try:
                unweighted_energy_score[idx] += mtx[seq[idx]][aa] * freq
            except KeyError:
                unweighted_energy_score[idx] += 0

    if mode == 'short':
        for idx in range(len(seq)):
            for idx2 in range(idx - wc, idx + wc + 1):
                if idx2 < 0 or idx2 >= len(seq):
                    weighted_energy_score[idx] += -1.26
                else:
                    weighted_energy_score[idx] += unweighted_energy_score[idx2]
            weighted_energy_score[idx] /= len(range(idx - wc, idx + wc + 1))
    else:
        weighted_energy_score = smooth(unweighted_energy_score, wc)

    glob_text = ""
    if mode == 'glob':
        gr = []
        in_gr = False
        beg, end = 0, 0
        for idx, val in enumerate(weighted_energy_score):
            if in_gr and val <= 0.3:
                gr.append({0: beg, 1: end})
                in_gr = False
            elif in_gr:
                end += 1
            if val > 0.3 and not in_gr:
                beg = idx
                end = idx
                in_gr = True
        if in_gr:
            gr.append({0: beg, 1: end})
        mgr = []
        k = 0
        kk = k + 1
        if gr:
            beg = gr[0][0]
            end = gr[0][1]
        nr = len(gr)
        while k < nr:
            if kk < nr and gr[kk][0] - end < 45:
                beg = gr[k][0]
                end = gr[kk][1]
                kk += 1
            elif end - beg + 1 < 35:
                k += 1
                if k < nr:
                    beg = gr[k][0]
                    end = gr[k][1]
            else:
                mgr.append({0: beg, 1: end})
                k = kk
                kk += 1
                if k < nr:
                    beg = gr[k][0]
                    end = gr[k][1]
        seq = seq.lower()
        nr = 0
        res = ""
        for i in mgr:
            res += seq[nr:i[0]] + seq[i[0]:i[1] + 1].upper()
            nr = i[1] + 1
        res += seq[nr:]
        res = " ".join([res[i:i + 10] for i in range(0, len(res), 10)])
        glob_text += "Number of globular domains: {}\n".format(len(mgr))
        for n, i in enumerate(mgr):
            glob_text += "          globular domain   {}.\t{}-{}\n".format(n + 1, i[0] + 1, i[1] + 1)
        glob_text += "\n".join(textwrap.wrap(res, 70))

    for idx, val in enumerate(weighted_energy_score):
        if val <= histo_min + 2 * histo_step:
            iupred_score[idx] = 1
        elif val >= histo_max - 2 * histo_step:
            iupred_score[idx] = 0
        else:
            iupred_score[idx] = histo[int((weighted_energy_score[idx] - histo_min) * (1 / histo_step))]
    return iupred_score, glob_text


def anchor2(folder, seq, iupred_scores):
    local_window_size = 41
    iupred_window_size = 30
    local_smoothing_window = 5
    par_a = 0.0013
    par_b = 0.26
    par_c = 0.43
    iupred_limit = par_c - (par_a / par_b)
    mtx = read_matrix(os.path.join(folder, 'anchor2_energy_matrix'))
    interface_comp = {}
    with open(os.path.join(folder, 'anchor2_interface_comp')) as _fn:
        for line in _fn:
            interface_comp[line.split()[1]] = float(line.split()[2])
    local_energy_score = [0] * len(seq)
    interface_energy_score = [0] * len(seq)
    energy_gain = [0] * len(seq)
    for idx in range(len(seq)):
        freq_dct = aa_freq(
            seq[max(0, idx - local_window_size):max(0, idx - 1)] + seq[idx + 2:idx + local_window_size + 1])
        for aa, freq in freq_dct.items():
            try:
                local_energy_score[idx] += mtx[seq[idx]][aa] * freq
            except KeyError:
                local_energy_score[idx] += 0
        for aa, freq in interface_comp.items():
            try:
                interface_energy_score[idx] += mtx[seq[idx]][aa] * freq
            except KeyError:
                interface_energy_score[idx] += 0
        energy_gain[idx] = local_energy_score[idx] - interface_energy_score[idx]
    iupred_scores = smooth(iupred_scores, iupred_window_size)
    energy_gain = smooth(smooth(energy_gain, local_smoothing_window), local_smoothing_window)
    anchor_score = [0] * len(seq)
    for idx in range(len(seq)):
        sign = 1
        if energy_gain[idx] < par_b and iupred_scores[idx] < par_c:
            sign = -1
        corr = 0
        if iupred_scores[idx] > iupred_limit and energy_gain[idx] < 0:
            corr = (par_a / (iupred_scores[idx] - par_c)) + par_b
        anchor_score[idx] = sign * (energy_gain[idx] + corr - par_b) * (iupred_scores[idx] - par_c)
        anchor_score[idx] = 1 / (1 + math.e ** (-22.97968 * (anchor_score[idx] - 0.0116)))
    return anchor_score
