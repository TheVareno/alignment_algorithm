class Alignment:
    def __init__(self, aligned_sequences, alignment_score):
        self.aligned_sequences = aligned_sequences
        self.alignment_score = alignment_score
    
    #def __str__(self):
     #   return self.aligned_sequences
                
    def get_aligned_sequences(self):
        return self.aligned_sequences
    
    def get_alignment_score(self): 
        return self.alignment_score
