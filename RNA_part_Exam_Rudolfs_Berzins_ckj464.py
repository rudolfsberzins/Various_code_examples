# Copyright 2012 Stefan E Seemann <seemann@rth.dk>
# Modified for assignment by Rudolfs Berzins <ckj464>
import numpy as np
import scipy.stats as stats

def pearsoncoff (wildtype_seq, mutant_seq):
    '''The function calculates Pearson's correlation coefficient by comparing Wildtype sequence with Mutant sequence.
     This is done by calculating their C-matrix, making a 1-D vector out of the matrix by suming rows and using this
     vector to calculate the correlation coefficient.'''

    def wildtype(seq_input):
        '''Function generates a Nussinov scoring matrix and a C-matrix, which then is used for calculation of a pi vector'''

        # input
        seq = seq_input
        seql = [i for i in seq.upper()]
        l = len(seq)

        # loop size
        h = 3

        # basepair scores
        scores = {'AU': 2, 'UA': 2, 'GU': 1, 'UG': 1, 'GC': 3, 'CG': 3}

        # initialize array
        m = [[0 for i in range(l)] for j in range(l)]

        # fill scoring matrix
        for j0 in range(h + 1, l):
            for i in range(0, l - j0):
                j = i + j0
                # rule 1) i,j paired
                if seql[i] + seql[j] in scores:
                    m[i][j] = m[i + 1][j - 1] + scores[seql[i] + seql[j]]
                # rule 2) i unpaired
                if m[i + 1][j] > m[i][j]:
                    m[i][j] = m[i + 1][j]
                # rule 3) j unpaired
                if m[i][j - 1] > m[i][j]:
                    m[i][j] = m[i][j - 1]
                # rule 4) bifurcation k
                for k in range(i + 1 + h, j - 1 - h):
                    if m[i][k] + m[k + 1][j] > m[i][j]:
                        m[i][j] = m[i][k] + m[k + 1][j]
                    # print (i,j,seql[i],seql[j],m[i][j])

        # backtracking
        str = ['.' for i in range(l)]
        stack = []
        stack.append([0, l - 1])
        c = np.zeros(shape=(l, l)) # initiate the C-matrix
        while len(stack) > 0:
            top = stack.pop(),
            i = top[0][0]
            j = top[0][1]
            if i >= j:
                continue
            elif m[i + 1][j] == m[i][j]:
                stack.append([i + 1, j])
            elif m[i][j - 1] == m[i][j]:
                stack.append([i, j - 1])
            elif seql[i] + seql[j] in scores and m[i + 1][j - 1] + scores[seql[i] + seql[j]] == m[i][j]:
                # record basepair i,j
                # print ("BP",i,j,":",seql[i],seql[j],"\n",)
                str[i] = "("
                str[j] = ")"
                stack.append([i + 1, j - 1])
                # does the backtracking and also fill in the C matrix if conditions for base pair are met
                if seql[i] + seql[j] in scores and m[i + 1][j - 1] + scores[seql[i] + seql[j]] == m[i][j] and seql[i]+seql[j] == 'AU':
                    c[i][j] = 2
                    c[j][i] = 2
                elif seql[i] + seql[j] in scores and m[i + 1][j - 1] + scores[seql[i] + seql[j]] == m[i][j] and seql[i]+seql[j] == 'UA':
                    c[i][j] = 2
                    c[j][i] = 2
                elif seql[i] + seql[j] in scores and m[i + 1][j - 1] + scores[seql[i] + seql[j]] == m[i][j] and seql[i]+seql[j] == 'GC':
                    c[i][j] = 3
                    c[j][i] = 3
                elif seql[i] + seql[j] in scores and m[i + 1][j - 1] + scores[seql[i] + seql[j]] == m[i][j] and seql[i]+seql[j] == 'CG':
                    c[i][j] = 3
                    c[j][i] = 3
                elif seql[i] + seql[j] in scores and m[i + 1][j - 1] + scores[seql[i] + seql[j]] == m[i][j] and seql[i]+seql[j] == 'GU':
                    c[i][j] = 1
                    c[j][i] = 1
                elif seql[i] + seql[j] in scores and m[i + 1][j - 1] + scores[seql[i] + seql[j]] == m[i][j] and seql[i]+seql[j] == 'AU':
                    c[i][j] = 1
                    c[j][i] = 1
            else:
                for k in range(i + 1 + h, j - 1 - h):
                    if m[i][k] + m[k + 1][j] == m[i][j]:
                        stack.append([k + 1, j])
                        stack.append([i, k])
                        break

        # calculate the pi vector
        for pi in range(len(c)):
            vec_pi_w = np.sum(c, axis=1)

        return vec_pi_w

    def mutant(seq_input):
        '''Function generates a Nussinov scoring matrix and a C-matrix, which then is used for calculation of a pi vector'''

        # input
        seq = seq_input
        seql = [i for i in seq.upper()]
        l = len(seq)

        # loop size
        h = 3

        # basepair scores
        scores = {'AU': 2, 'UA': 2, 'GU': 1, 'UG': 1, 'GC': 3, 'CG': 3}

        # initialize array
        m = [[0 for i in range(l)] for j in range(l)]

        # fill scoring matrix
        for j0 in range(h + 1, l):
            for i in range(0, l - j0):
                j = i + j0
                # rule 1) i,j paired
                if seql[i] + seql[j] in scores:
                    m[i][j] = m[i + 1][j - 1] + scores[seql[i] + seql[j]]
                # rule 2) i unpaired
                if m[i + 1][j] > m[i][j]:
                    m[i][j] = m[i + 1][j]
                # rule 3) j unpaired
                if m[i][j - 1] > m[i][j]:
                    m[i][j] = m[i][j - 1]
                # rule 4) bifurcation k
                for k in range(i + 1 + h, j - 1 - h):
                    if m[i][k] + m[k + 1][j] > m[i][j]:
                        m[i][j] = m[i][k] + m[k + 1][j]
                    # print (i,j,seql[i],seql[j],m[i][j])

        # backtracking
        str = ['.' for i in range(l)]
        stack = []
        stack.append([0, l - 1])
        c = np.zeros(shape=(l, l)) # initiate the C-matrix
        while len(stack) > 0:
            top = stack.pop(),
            i = top[0][0]
            j = top[0][1]
            if i >= j:
                continue
            elif m[i + 1][j] == m[i][j]:
                stack.append([i + 1, j])
            elif m[i][j - 1] == m[i][j]:
                stack.append([i, j - 1])
            elif seql[i] + seql[j] in scores and m[i + 1][j - 1] + scores[seql[i] + seql[j]] == m[i][j]:
                # record basepair i,j
                # print ("BP",i,j,":",seql[i],seql[j],"\n",)
                str[i] = "("
                str[j] = ")"
                stack.append([i + 1, j - 1])
                # does the backtracking and also fill in the C matrix if conditions for base pair are met
                if seql[i] + seql[j] in scores and m[i + 1][j - 1] + scores[seql[i] + seql[j]] == m[i][j] and seql[i]+seql[j] == 'AU':
                    c[i][j] = 2
                    c[j][i] = 2
                elif seql[i] + seql[j] in scores and m[i + 1][j - 1] + scores[seql[i] + seql[j]] == m[i][j] and seql[i]+seql[j] == 'UA':
                    c[i][j] = 2
                    c[j][i] = 2
                elif seql[i] + seql[j] in scores and m[i + 1][j - 1] + scores[seql[i] + seql[j]] == m[i][j] and seql[i]+seql[j] == 'GC':
                    c[i][j] = 3
                    c[j][i] = 3
                elif seql[i] + seql[j] in scores and m[i + 1][j - 1] + scores[seql[i] + seql[j]] == m[i][j] and seql[i]+seql[j] == 'CG':
                    c[i][j] = 3
                    c[j][i] = 3
                elif seql[i] + seql[j] in scores and m[i + 1][j - 1] + scores[seql[i] + seql[j]] == m[i][j] and seql[i]+seql[j] == 'GU':
                    c[i][j] = 1
                    c[j][i] = 1
                elif seql[i] + seql[j] in scores and m[i + 1][j - 1] + scores[seql[i] + seql[j]] == m[i][j] and seql[i]+seql[j] == 'AU':
                    c[i][j] = 1
                    c[j][i] = 1
            else:
                for k in range(i + 1 + h, j - 1 - h):
                    if m[i][k] + m[k + 1][j] == m[i][j]:
                        stack.append([k + 1, j])
                        stack.append([i, k])
                        break

        # calculate the pi vector
        for pi in range(len(c)):
            vec_pi_m = np.sum(c, axis=1)
        print (seq,l,"\n",''.join(str),"\n","Score:",m[0][l-1])
        return vec_pi_m

    # Generate vectors for Wildtype and Mutant sequence
    wildtype_vector = wildtype(wildtype_seq)
    mutant_vector = mutant(mutant_seq)

    # Calculate Pearson's correlation coefficient
    coff = stats.pearsonr(wildtype_vector,mutant_vector)

    print(coff)

# Used sequences
Wild = 'ggacgaaucucuggagagacucccucucgcuuuaaauagcguagaggaaaacgagcaccgaaggagcaaauccgcuacuauagcggauaaucucucagguaaaaggacagagacaagcga'
Mutant1 = 'ggacgaaucucuggagagacucccucucgcuuuaaauagcguagaggaaaacgagcaccgaaggagcaaauccgcuacuauagcggCuaaucucucagguaaaaggacagagacaagcga'
Mutant2 = 'ggacgaaucucuggagagacucUcucucgcuuuaaauagcguagaggaaaacgagcaccgaaggagcaaauccgcuacuauagcggauaaucucucagguaaaaggacagagacaagcga'
Mutant3 = 'ggacgaaucucuggagagacucccucucgcuuuaaauagcguagaggaaaacgagcaccgaaggagcaaauccgcuCcuauagcggauaaucucucagguaaaaggacagagacaagcga'
Mutant4 = 'ggacgaaucucuggagagacucccucucgcuuuaaauagcguagaggaaaacgagcaccgaaggagcaaauccgcuGcuauagcggauaaucucucagguaaaaggacagagacaagcga'

# Initiate calculation
pearsoncoff(Wild, Mutant1)
pearsoncoff(Wild, Mutant2)
pearsoncoff(Wild, Mutant3)
pearsoncoff(Wild, Mutant4)
pearsoncoff(Wild, Wild)




