# -*- coding: utf-8 -*-
"""
Created on Fri Mar 3 2023 16:24:15

@author: Fangrui Hu
"""
import time
import re
import math
import os
import sys
import numpy as np
import pandas as pd
import struct
from prody import *
import argparse
parser = argparse.ArgumentParser(description="Input parameters")
# Main parameters
parser.add_argument("--PDBID", type=str, help="Four-character pdbid (no .pdb required)", required=True)
parser.add_argument("--CHAIN", type=str, default=False, help="Protein chain (optional)")
sys.setrecursionlimit(100000)


class AllosESutils:
    def __init__(self, path, pdbid, chain=None):
        self.path = path
        self.pdbid = pdbid.upper()
        self.chain = chain.upper()

    def fpocket_run(self):
        os.system('fpocket ' + '-f ' + self.path + self.pdbid.upper() + '.pdb')

    def pdbchain(self):
        pdbname = self.pdbid
        chain = self.chain.upper()
        with open(self.path + pdbname.lower() + '.pdb', 'r') as f:
            pdb = f.readlines()
            for i in pdb:
                if i[:4] == 'ATOM' and (i[20:22].strip() == chain) and i[16] != 'B' and i[16] != 'C':
                    with open(self.path + pdbname.upper() + '_' + chain + '.pdb', 'a') as Chain:
                        Chain.write(i)
        return pdbname.upper() + '_' + chain

    def remove(self):
        os.system('chmod 777 ' + self.path + self.pdbid.upper() + '.pdb')
        os.system('rm' + ' ' + '-f' + ' ' + self.path + self.pdbid.upper() + '.pdb')
    # INFO
    def ext_info(self):
        # Physicochemical properties
        fpocket_info = ['Score :', 'Druggability Score :', 'Number of Alpha Spheres :', 'Total SASA :', 'Polar SASA :',
                        'Apolar SASA :', 'Volume :',
                        'Mean local hydrophobic density :', 'Mean alpha sphere radius :',
                        'Mean alp. sph. solvent access :', 'Apolar alpha sphere proportion :',
                        'Hydrophobicity score:', 'Volume score:', 'Polarity score:', 'Charge score :',
                        'Proportion of polar atoms:', 'Alpha sphere density :',
                        'Cent. of mass - Alpha Sphere max dist:', 'Flexibility :']
        # Read arguments
        input_file = self.path + self.pdbid
        input_file += '.pdb'  # Read PDB file
        proPDB = parsePDB(input_file)
        calphas = proPDB.select('calpha')
        # Extract fpocket_info
        pockets_info = self._extinfo(fpocket_info)
        # Extract pocket residues into a dictionary
        pockres, pockres_name = self._extpockets( calphas)
        clean_string = ['fpocket_' + re.sub('[:.\-\s+]', '', info_id) for info_id in fpocket_info] + ['residues']
        info = [[str(pockets_info[pockidx][info_id]) for info_id in fpocket_info] + [
            ','.join(pockres_name[pockidx])] for pockidx in range(len(pockres))]
        info = pd.DataFrame(info, columns=clean_string).drop(columns=['fpocket_DruggabilityScore',
                                                                      'fpocket_NumberofAlphaSpheres',
                                                                      'fpocket_TotalSASA', 'fpocket_PolarSASA',
                                                                      'fpocket_ApolarSASA', 'fpocket_Volume',
                                                                      'fpocket_Apolaralphasphereproportion',
                                                                      'fpocket_Hydrophobicityscore',
                                                                      'fpocket_Chargescore',
                                                                      'fpocket_Alphaspheredensity',
                                                                      'fpocket_CentofmassAlphaSpheremaxdist', ])
        time.sleep(0.2)
        return info
    #DSSP
    def dssp(self):
        dssp_format = '5s5s2s2s3s4s2s15s12s11s11s11s8s6s6s6s6s7s7s7s1s'
        SSEfeature = []
        for pocketnum in range(len(self.ext_pocket())):
            pocket = self.ext_pocket()[pocketnum]
            dssp = open(self.path + self.pdbid + '.dssp', 'r')  # Open the dssp file
            line = dssp.readlines()
            title_index = 0
            H_num = 0
            B_num = 0
            E_num = 0
            G_num = 0
            I_num = 0
            T_num = 0
            S_num = 0
            M_num = 0
            for i in range(len(line)):
                if line[i][0:3].strip() == '#':
                    title_index = i  # Determine the row index of the dssp header
            for i in range(title_index + 1, len(line)):
                col = struct.unpack(dssp_format, line[i].encode())
                type = col[4].strip().decode("utf-8")
                resname = col[3].strip().decode("utf-8")
                resindex = col[1].strip().decode("utf-8")
                reschain = col[2].strip().decode("utf-8")
                if (resindex + ':' + reschain) in pocket:
                    if resname != '!':
                        if type == 'H':
                            H_num = H_num + 1
                        if type == 'B':
                            B_num = B_num + 1
                        if type == 'E':
                            E_num = E_num + 1
                        if type == 'G':
                            G_num = G_num + 1
                        if type == 'I':
                            I_num = I_num + 1
                        if type == 'T':
                            T_num = T_num + 1
                        if type == 'S':
                            S_num = S_num + 1
                        if type == '':
                            M_num = M_num + 1
            residue_num = H_num + B_num + E_num + G_num + I_num + T_num + S_num + M_num
            H_per = (H_num / residue_num) * 100

            SSEfeature.append(
                [self.pdbid + '_pocket' + str(pocketnum + 1), str(H_per), str(E_num), str(G_num), str(M_num)])
        SSEfeatures = np.array(SSEfeature)
        SSE = pd.DataFrame(SSEfeatures,
                           columns=['pocket_name', 'pocket_H_per', 'pocket_E_num', 'pocket_G_num', 'pocket_M_num'])
        time.sleep(0.2)
        return SSE

    def ext_pocket(self):
        # Read arguments
        input_prefix = self.path + self.pdbid
        pdb_file = input_prefix + '.pdb'
        # Read PDB file
        pro = parsePDB(pdb_file)
        calphas = pro.select('calpha')
        pockets_res, pockets_res_name = self._extpockets(calphas)
        return pockets_res_name


    # [Allopred]
    def _extinfo(self, fpocket_info):
        pockinfo = {}
        pockfile = open(self.path + self.pdbid + '_out/' + self.pdbid + '_info.txt')
        pocklines = pockfile.readlines()
        pockfile.close()
        # Create dictionaries
        pockcounter = 0
        for line in pocklines:
            if line[:6] == 'Pocket':
                pockinfo[pockcounter] = {}
                pockcounter += 1
        # Write to dictionary
        for fp_info in fpocket_info:
            count = 0
            for line in pocklines:
                if line.lstrip()[:len(fp_info)] == fp_info:
                    split = re.split(r'\s+', line.rstrip())
                    pockinfo[count][fp_info] = float(split[-1])
                    count += 1
        return pockinfo

    def _extpockets(self, selection):
        pockcount = 0  # Pocket counter
        pockfile = open(self.path + self.pdbid + '_out/' + self.pdbid + '_info.txt')
        for line in pockfile:
            if line[:6] == 'Pocket':
                pockcount += 1
        pockfile.close()
        pockres = {}  # Pocket Index
        pockres_name = {}  # Residue number and chain
        for pockidx in range(pockcount):
            pockres[pockidx] = []
            pockres_name[pockidx] = []
            in_file = open(self.path + self.pdbid + '_out/pockets/pocket' + str(pockidx) + '_atm.pdb')
            for line in in_file:
                if line[:4] == 'ATOM':
                    # Add residue name
                    resname = line[22:26].lstrip() + ':' + line[21]
                    if resname not in pockres_name[pockidx]:
                        pockres_name[pockidx].append(resname)
                    # Add residue ID
                    selector = 'chain ' + line[21] + ' and resnum ' + line[22:26].lstrip()
                    try:
                        index = selection.select(selector).getResindices()[0]
                        if index not in pockres[pockidx]:
                            pockres[pockidx].append(index)
                    except TypeError:
                        pass
                    except AttributeError:
                        pass
                    except atomic.select.SelectionError:
                        pass
            in_file.close()
        return pockres, pockres_name

    # SNB-PSSM
    # Visualization
    # Extract sequence annotation program
    def sequenceComment(self, line):
        t = 0
        comment = ''
        for j in line:
            if j == '"':
                t += 1
            if t == 1:
                comment = comment + j
            else:
                if t == 2:
                    return comment

    # Extraction sequence encoding procedure
    def sequenceCode(self, sequence):
        change_str = ''
        for i in sequence:
            i.strip()
            change_str = change_str + i
        t = 0
        outSequence = ''
        for j in change_str:
            if j == '"':
                t += 1
            if t == 1:
                outSequence = outSequence + j
            else:
                if t == 2:
                    return outSequence[1:]

    # Extraction of PSSM matrix
    def pssm_matrix(self, d, x, pssmBox):
        num = d[x + 3:x + 23]
        num1 = d[x + 1]
        outLine = ''
        count = 0
        for j in num:
            count = count + 1
            if count == 19:
                line = num1.replace(',\n', '')
                outLine = outLine + line
            else:
                line = j.replace(',\n', '')
                outLine = outLine + line
        pssmBox.append(outLine)
        if x + 28 <= len(d) - 30:
            call = self.pssm_matrix(d, x + 28, pssmBox)
        return pssmBox

    # Single row PSSM matrix rearrangement
    def singlePSSM(self, pssmBox, rowName):
        midBox = []
        pssmLine = ''
        pssmRow = pssmBox[rowName].strip()
        midPSSM = pssmRow.replace('        ', ' ')
        midBox = midPSSM.split(' ')
        for j in midBox:
            if len(str(j)) == 1:
                newNumber = ' ' + j
                pssmLine = pssmLine + '  ' + newNumber
            else:
                newNumber = j
                pssmLine = pssmLine + '  ' + newNumber
        return pssmLine

    # Auxiliary marking procedures
    def help_Pro(self, d):
        teq = ['1', ]
        for i in range(len(d)):
            if 'seq-data' in d[i]:
                teq.clear()
            if len(teq) == 0:
                if '}' in d[i]:
                    d[i] = '}12345'
                    return d
                else:
                    continue

    def pssm_pro(self, formerName):
        with open(formerName, 'r') as u:
            d = u.readlines()
            u.close()
        print('Read successfully, PSSM matrix being extracted...')
        d = self.help_Pro(d)
        seq = ['1', ]
        sequence = ''
        yourwant = []
        pssmBox = []
        time.sleep(0.2)
        for line in range(len(d)):
            i = d[line]
            i.strip()
            if 'seq-data' in i:
                seq.clear()
                seq.append(i)
                yourwant.extend(seq)
                seq.clear()
            if '}12345' in i:
                sequence = yourwant[1:]
                seq = ['1', ]
                sequence = self.sequenceCode(sequence)
                sequence = sequence.replace('\n', '')
            else:
                if len(seq) == 0:
                    yourwant.append(i)
            if 'finalData' in i:
                pssmBox = self.pssm_matrix(d, line + 2, pssmBox)  # Extraction of PSSM matrix
        print('Extraction successful, exporting file now...')
        return pssmBox, sequence

    def PSSM_matrix(self, sequence, pssmBox):
        Result = ''
        for rowName in range(len(sequence)):
            pssmLine = self.singlePSSM(pssmBox, rowName)  # Rearrange single row PSSM matrix
            trueName = rowName + 1
            if len(str(trueName)) == 1:
                rowLine = '  ' + str(trueName)
            else:
                if len(str(trueName)) == 2:
                    rowLine = ' ' + str(trueName)
                else:
                    rowLine = str(trueName)
            outRow = rowLine + ' ' + sequence[rowName] + '  ' + pssmLine
            Result = Result + '\n' + outRow

        title = ['NUM', 'RES', 'C', 'D', 'E ', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V',
                 'W',
                 'A ', 'Y ']
        dataxlsx = []
        for i in range(len(Result.split('\n')) - 1):
            for j in range(len(Result.split('\n')[i + 1].split())):
                dataxlsx.append(Result.split('\n')[i + 1].split()[j])
        dataxlsx = np.array(dataxlsx).reshape(len(Result.split('\n')) - 1, 22)

        dataxlsx = pd.DataFrame(dataxlsx, columns=title)
        for c in dataxlsx.columns:
            if c != 'RES':
                try:
                    dataxlsx[c] = dataxlsx[c].astype(np.int64)
                except:
                    print('error')
        return dataxlsx

    # Calculating SNB-PSSM
    def SNB_PSSM(self, data):
        posall = []
        residx = []
        chain = []
        with open(self.path + self.pdbid + '.pdb', 'r') as fid1:
            fid1 = fid1.readlines()
        for pl in fid1:
            pl = pl.strip('\n')
            if len(pl) < 3:
                break
            elif pl[0:4] == 'ATOM' and pl[11:15].strip() == 'CA' and pl[16] != 'B' and pl[16] != 'C':
                residx.append(pl[22:26].strip())
                chain.append(pl[20:22].strip())
                posall.append([float(pl[31:38].strip()), float(pl[39:46].strip()), float(pl[47:54].strip())])

        posall = np.array(posall)
        # CA atoms within 7.5 Å are calculated and these residues are summed and averaged
        row = []  # Start Node
        col = []  # End Node
        for i in range(len(posall)):
            connected_index = []
            for j in range(len(posall)):
                distance = np.sqrt(np.sum(np.square(posall[i, :] - posall[j, :])))
                if distance < 7.5:
                    connected_index.append(j)
            if connected_index:
                col.extend([i] * len(connected_index))
                row.extend(connected_index)
        pssm = data
        res = pssm['RES'].to_numpy()
        pssm = pssm.drop(columns=['NUM', 'RES'])
        D = pssm.iloc[row, :]
        k = pd.DataFrame(col).diff()
        k = k.loc[k.values != 0].index.values
        PSSM = [np.mean(D.iloc[k[i - 1]:k[i], :]).tolist() for i in range(1, len(k))]
        PSSM.append(np.mean(D.iloc[k[-1]:len(col), :]).tolist())
        save = pd.DataFrame(np.c_[residx, chain, res, PSSM],
                            columns=['residx', 'chain', 'resname'] + pssm.columns.tolist())

        return save

    # Extraction Pocket SNB-PSSM
    def pocket_SNB(self):
        pssmBox, sequence = self.pssm_pro(self.path + self.pdbid + '.asn')
        datacsv = self.PSSM_matrix(sequence, pssmBox)  # PSSM
        PSSM = self.SNB_PSSM(datacsv)  # SNB-PSSM
        PSSMmat = []
        pocketre = self.ext_pocket()
        # The residue number corresponds to the net transfer entropy put into the dictionary
        Mappres = {str(PSSM.iloc[q, 0]) + ':' + str(PSSM.iloc[q, 1]): PSSM.iloc[q, 3:].tolist() for q in
                   range(len(PSSM))}
        # Store the data, extract the pocket residues, and take the mean value
        for pocket_key in pocketre:
            pocketid = self.pdbid + '_pocket' + str(pocket_key + 1)  # Pocket Index
            resid = pocketre[pocket_key]  # Pocket residual base number
            if pocketid[:4] == self.pdbid[:4]:
                val = [Mappres[res] for res in resid]
                PSSMmat.append(np.array(val).mean(axis=0).tolist())
        pocket_SNB = pd.DataFrame(np.array(PSSMmat),
                                  columns=['C', 'D', 'E ', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S',
                                           'T',
                                           'V', 'W', 'A', 'Y']).drop(
            columns=['E ', 'F', 'G', 'I', 'L', 'N', 'P', 'Q', 'T', 'V', 'Y'])
        time.sleep(0.2)
        return pocket_SNB

    def read_pdb(self):
        pdb = open(self.path + self.pdbid + '.pdb', 'r')
        residues_index = []
        residues_name = []
        residues_coordinate = []
        residues_B_factor = []
        residues_chain = []
        for atom in pdb:
            if len(atom) <= 3:
                continue
            elif atom[0:3] == '':
                break
            elif ((atom[0:4] == 'ATOM')
                  and atom[13:15] == 'CA' and atom[16] != 'B' and atom[16] != 'C'):
                coordinate_x = float(atom[30:38])
                coordinate_y = float(atom[38:46])
                coordinate_z = float(atom[46:54])
                residues_coordinate.append([coordinate_x, coordinate_y, coordinate_z])
                residues_B_factor.append(float(atom[60:68]))
                residues_index.append(int(atom[22:26]))
                residues_name.append(atom[17:20])
                residues_chain.append(atom[21:22])
                continue
        pdb.close()
        N_residues = len(residues_index)

        return N_residues, residues_index, residues_name, residues_coordinate, residues_chain

    def distance(self, coordinate_matrix):
        size, _ = np.shape(coordinate_matrix)
        dis = np.zeros([size, size])
        for i in range(size):
            for j in range(size):
                if j == i:
                    continue
                else:
                    dis[i, j] = math.sqrt((coordinate_matrix[i, 0] - coordinate_matrix[j, 0]) ** 2
                                          + (coordinate_matrix[i, 1] - coordinate_matrix[j, 1]) ** 2
                                          + (coordinate_matrix[i, 2] - coordinate_matrix[j, 2]) ** 2)
        return dis

    def kirchhoff(self, coordinate_matrix, cutoff):
        size, _ = np.shape(coordinate_matrix)
        Kirchhoff_matrix = np.zeros([size, size])
        dis = self.distance(coordinate_matrix)
        for i in range(size):
            for j in range(size):
                if j == i:
                    continue
                elif j != i:  # Non-diagonal
                    if dis[i, j] <= cutoff:
                        Kirchhoff_matrix[i, j] = -1
                    else:
                        Kirchhoff_matrix[i, j] = 0
            Kirchhoff_matrix[i, i] = -1 * sum(Kirchhoff_matrix[i, :])
        return Kirchhoff_matrix

    def GNM(self, cutoff):
        [N_P, _, _, residues_coordinate, _] = self.read_pdb()
        N = N_P
        coordinate = np.array(residues_coordinate)
        Kirchhoff = self.kirchhoff(coordinate, cutoff)
        [Vectors, Values, VectorsT1] = np.linalg.svd(Kirchhoff)
        sorted_indices = np.argsort(Values)
        Values = Values[sorted_indices[:: 1]]
        Vectors = Vectors[:, sorted_indices[:: 1]]
        InvKirchhoff = (Vectors) * (np.linalg.pinv(np.diag(Values))) * (Vectors.T)
        CellAij = {}
        for k in range(0, N):
            if (1 / Values[k]) < 1000:
                CellAij[k] = (Vectors[:, k] * (np.array([Vectors[:, k]]).T) / Values[k])
            else:
                CellAij[k] = np.zeros([N, N])
        return InvKirchhoff, CellAij, N, Values

    def Transfer_entropy(self, cutoff, tau):
        InvKirchhoff, CellAij, N, eig_Value = self.GNM(cutoff)
        TE = np.ones((N, N), dtype=np.complex_)
        for i in range(N):
            for j in range(N):
                aEk = [CellAij[k][j][j] for k in range(0, N)]
                bEk = [CellAij[k][i][j] for k in range(0, N)]
                cEk = [CellAij[k][j][j] for k in range(0, N)]
                dEk = [CellAij[k][i][j] for k in range(0, N)]
                eEk = [CellAij[k][i][i] for k in range(0, N)]
                aEk = aEk * np.exp(-eig_Value * tau)
                bEk = bEk * np.exp(-eig_Value * tau)
                a = np.sum(cEk) ** 2 - np.sum(aEk) ** 2
                b = (np.sum(eEk) * np.sum(cEk) ** 2)
                c = 2 * (np.sum(dEk)) * np.sum(aEk) * np.sum(bEk)
                d = -(((np.sum(bEk) ** 2) + (np.sum(dEk) ** 2)) * (np.sum(cEk))) - ((np.sum(aEk) ** 2) * np.sum(eEk))
                f = np.sum(cEk)
                g = (np.sum(eEk) * np.sum(cEk)) - (np.sum(dEk) ** 2)
                if i == j:
                    TE[i][j] = 0
                else:
                    TE[i][j] = 0.5 * np.log(a) - 0.5 * np.log(b + c + d) - 0.5 * np.log(f) + 0.5 * np.log(g)
        TE[TE < 0] = 0
        netTE = TE - TE.T
        Difference = np.real((netTE).sum(axis=1))
        norm_difference = Difference / np.max(np.abs(Difference))
        return norm_difference

    # NTE
    def pocket_NTE(self):
        cutoff = 7
        tau = 5
        NTE = self.Transfer_entropy( cutoff, tau)
        pocketre = self.ext_pocket()
        # All residue numbers of the pdb file
        _, residues_index, _, _, chain = self.read_pdb()
        # The residue number corresponds to the net transfer entropy put into the dictionary
        Mappres = {str(residues_index[q]) + ':' + chain[q]: NTE[q] for q in range(len(NTE))}
        NTEmat = []
        # Store the data, extract the pocket residues, and take the mean value
        for pocket_key in pocketre:
            pocketid = self.pdbid + '_pocket' + str(pocket_key + 1)  # Pocket Index
            resid = pocketre[pocket_key]  # Pocket residual base number
            if pocketid[:4] == self.pdbid[:4]:
                val = [float(Mappres[res]) for res in resid]
                val.sort(key=lambda x: float(x), reverse=True)
                NTEmat.append(np.average(val))
        pocket_NTE = pd.DataFrame(np.array(NTEmat), columns=['ave_NTE'])
        time.sleep(0.2)
        return pocket_NTE

    def original_features(self):
        print('******Detection pockets...******')
        self.fpocket_run()
        DSSP = self.dssp()
        NTE = self.pocket_NTE()
        PSSM_SNB = self.pocket_SNB()
        INFO = self.ext_info()
        summary = pd.concat([DSSP, NTE, PSSM_SNB, INFO], axis=1)
        print('********************************************')
        print('******Program starts running...******')
        print('******Start extracting features...******')
        print('******Feature integration completed!******')
        return summary
