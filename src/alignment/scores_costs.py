
from alignment.sequence import Sequence

class ScoresCosts:
    def __init__(self, alignment_type, sequence_type):
        self.alignment_type = alignment_type
        self.sequence_type = sequence_type
        self.score_cost = None    
    
    def unit_costs(self, a, b):
        if a == b and self.alignment_type == 'global':
            return 0
        elif a == b and self.alignment_type == 'local':
            return 1
        elif a != b and self.alignment_type == 'global':
            return 1
        elif a != b and self.alignment_type == 'local':
            return 0 
    
    # base on blosum matrix costs     
    def blosum(a, b):
        # define blosum as a dictionary : keys are tuple from two amino acids and values are corresponding score
        blosum_matrix = {
            ('C', 'C') : 9, ('C', 'S'): -1, ('C', 'T') : -1, ('C', 'A') : 0, ('C', 'G') : -3 , ('C', 'P') : -3, 
            ('C', 'D') : -3, ('C', 'E') : -4, ('C', 'Q') : -3, ('C', 'N') : -3, ('C', 'H') : -3, ('C', 'R') : -3, ('C', 'K') : -3, 
            ('C', 'M') : -1, ('C', 'I') : -1, ('C', 'L') : -1, ('C', 'V'): -1, ('C', 'W') : -2, ('C', 'Y') : -2, ('C', 'F') : -2,
            ('S', 'S') : 4, ('S', 'T'): 1, ('S', 'A') : 1, ('S', 'G') : 0, ('S', 'P') : -1, ('S','D') : -1, ('S', 'E') : 0, ('S', 'Q'): 0,
            ('S', 'N') : 1, ('S', 'H') : -1, ('S', 'R'): -1, ('S', 'R') : -1, ('S', 'K') : 0,
            ('S', 'M') : -1, ('S', 'I') : -2, ('S', 'L') : -2, ('S', 'V') : -2, ('S', 'W') : -3,
            ('S', 'W') : -2, ('S', 'F') : -2, ('T', 'T') : 5, ('T', 'A'): 0, ('T', 'G') : -2, ('T', 'P') : -1, ('T', 'D') : -1, ('T', 'E') : -1, ('T', 'Q'): -1, ('T', 'N'): 0,
            ('T', 'H') : -2, ('T', 'R') : -1, ('T', 'K') : -1, ('T', 'M') : -1, ('T', 'I') : -1, ('T', 'L') : -1, ('T', 'V') : 0, ('T', 'W') : -2, ('T', 'Y') : -2, ('T', 'F') : -2,
            ('A', 'A') : 4, ('A', 'G') : 0, ('A', 'P') : -1, ('A', 'D') : -2, ('A', 'E') : -1, ('A', 'Q') : -1, 
            ('A', 'N'): -2, ('A', 'H'): -2, ('A', 'R') : -1, ('A', 'K') : -1, ('A', 'M') : -1, ('A', 'I'): -1, 
            ('A', 'L') : -1, ('A', 'V'): 0, ('A', 'W') : -3, ('A', 'Y') : -2, ('A', 'F') : -2,
            ('G', 'G') : 6, ('G', 'P') : -2, ('P', 'P'): 7, ('G', 'D') : -1, ('G','E') : -2, ('G','Q'): -2, ('G', 'N') : 0, 
            ('P', 'D') : -1, ('P', 'E') : -1, ('P', 'Q') : -1, ('P', 'N') : -2, ('G', 'H') : -2, ('G', 'R') : -2, ('G', 'K'): -2,
            ('P', 'H') : -2, ('P', 'R') : -2, ('P', 'K') : -1 
            
        }
        
       

    
    def set_score_cost(self): 
        if self.sequence_type == 'nt':
            self.unit_costs()
            
        elif self.sequence_type == 'aa':
            self.blosum()        