
from src.alignment.sequence import Sequence as seq
from src.alignment.alignment import Alignment 
from src.alignment.scores_costs import ScoresCosts  
 
class Algorithms(ScoresCosts): 
    def __init__(self, alignment_type, sequence_type):
        super().__init__(alignment_type, sequence_type)
        
    # smith waterman algorithm 
    def local_alignment(self, sequence1, sequence2):
        seq_type = self.sequence_type
        # initilze the DP Matrix  
        sequence1 = str(sequence1)
        sequence2 = str(sequence2)
        
        # create a ziro-filled matrix 
        score = [[0 for j in range(len(sequence2) + 1)] for i in range(len(sequence1) + 1)]
        
        # initilize the traceback matrix 
        traceback = [[0 for j in range(len(sequence2) + 1 )] for i in range(len(sequence1) + 1)]
        
        max_score = 0
        max_i = 0
        max_j = 0
        # start filling matrix and aligning the sequences based on sequence type 
        if seq_type == 'nt':
            for i in range(len(sequence1)):
                for j in range(len(sequence2)):
                    # To avoid IndexError put zeros index by i and j manually 
                    match_score = 0
                    insert_score = 0
                    delete_score = 0
                    if i == 0 and j == 0 :
                        match_score = 0
                    elif i == 0 and j != 0 :
                        insert_score = score[0][j]  - 1 
                    elif i != 0 and j == 0 :
                        delete_score = score[i][0]  - 1 
                    else:
                        match_score = score[i-1][j-1] + self.unit_costs(sequence1[i-1], sequence2[j-1])
                        delete_score = score[i-1][j]  - 1 
                        insert_score = score[i][j-1]  - 1

                    # now fill the current cell with the maximum of three possible calculation  
                    score[i][j] = max(match_score, delete_score, insert_score)
                    
                    # start filling traceback matrix in order to perform backtracking  
                    if score[i][j] == match_score:
                        traceback[i][j] = 1
                    elif score[i][j] == delete_score:
                        traceback[i][j] = 2
                    else:
                        traceback[i][j] = 3 

                if score[i][j] > max_score:
                    max_score = score[i][j]
                    max_i, max_j = i, j
                    
            # determine the optimal aligniment bw seq1 and seq2 : backtracking based on trace matrix 
                ui , vj = '', ''
                i, j = max_i , max_j

                # start going back from the maximum number to the start point 
                while score[i][j] > 0:
                    if traceback[i][j] == 1 : 
                        ui = sequence1[i-1] + ui
                        vj = sequence2[j-1] + vj
                    # a diagonal backing movement in DP Matrix :
                        i, j = i-1, j-1 
                    # upside movement : insertion
                    elif traceback[i][j] == 3: 
                        ui = sequence1[i-1] + ui
                        vj = '-' + vj
                        i -= 1
                    # leftside movement : deletion
                    elif traceback[i][j] == 2: 
                        ui = '-' + ui
                        vj = sequence2[j-1] + vj
                        j -= 1
        
        # now calculation of similarity based on blosum matrix
        elif seq_type == 'aa':
            pass
        
        alined_sequences = f"{ui}+\n+{vj}"
        global alignment 
        alignment = Alignment(alined_sequences ,max_score)
        return alignment
    
    # needleman wunsch algorithm              
    def global_alignment(self, sequence1, sequence2):
        

        # when the input sequences are dna 
        if self.alignment_type == 'dna':
            
            # define the DP matrix for global alignment
            cost = [[0 for j in range(len(sequence2))] for i in range(len(sequence1))] 

            # define the traceback matrix in order to determine the global optimal alignment
            traceback = [[0 for j in range(len(sequence2))] for i in range(len(sequence1))]
            
            # fill the DP matrix 
            for i in range(len(sequence1)):
                for j in range(len(sequence2)):
                    if i == 0 and j == 0:
                        match_cost = 0
                    elif i == 0 and j != 0 :
                        insert_cost = cost[0][j-1] + 1
                    elif i != 0 and j == 0 : 
                        delete_cost = cost[i-1][0] + 1  
                    else : 
                        match_cost = cost[i-1][j-1] + self.unit_costs(sequence1[i-1], sequence2[j-1])
                        insert_cost = cost[i-1][j] + 1
                        delete_cost = cost[i][j-1] + 1
                    # now fill the current cell of the matrix with minimum value of last cells 
                    cost[i][j] = min(match_cost, delete_cost, insert_cost)    
                    
                    # parallel we fill the trace matrix 
                    if cost[i][j] == match_cost:
                        traceback[i][j] == 1
                    elif cost[i][j] == delete_cost : 
                        traceback[i][j] == 2
                    elif cost[i][j] == insert_cost:
                        traceback[i][j] == 3
            
            u, v = '', ''
            k = len(sequence1)
            l = len(sequence2)
            end_cost = cost[k][l]
            # determine the alignment 
            while k > 0 and l > 0:
                if traceback[k][l] == 1:
                    u = u + sequence1[k]
                    v = v + sequence2[l]
                    k = k - 1
                    l = l - 1
                if traceback[k][l] == 2:
                    u = '-' + u 
                    v = sequence2[l] + v
                    k = k - 1
                if traceback[k][l] == 3:
                    u = sequence1[k] + u
                    v = '-' + v 
                    l = l - 1
            
            while l > 0:
                u = sequence1[l] + u
                v = '-' + v
            
            while k > 0:
                u = '-' + u
                v = sequence1[l] + v
            
            aligned_sequences = ''
            global alignment 
            alignment = Alignment(aligned_sequences, end_cost)             
            return alignment
        # when the input sequences are protein 
        if self.alignment_type == 'protein' :
            pass
        

alg = Algorithms('local', 'nt')
alg.alignment_type
alg.sequence_type
t = seq('AACCGGTT', 'nt')
t.sequence
t.type
z = seq('GTT', 'nt')


align = alg.local_alignment(t, z)
type(align)
align.aligned_sequences
align.alignment_score
