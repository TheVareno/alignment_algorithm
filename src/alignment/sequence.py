#/opt/anaconda3/bin/python

class Sequence: 
    def __init__(self, sequence_data, type):
        self.sequence = sequence_data # chars of the sequence
        self.type = type #nucleotide or aminoacids
        
    # convert the sequence object to string 
    def __str__(self):
        return self.sequence
    
    # return the length of the sequence object 
    def get_length(self): 
       return len(self.sequence)
    
    # return the squence object
    def get_sequence(self):
        return self.sequence
    
    def validate_sequence(self):
        pass
    
    def calculate_GC_content(self):
        pass
    
    def get_reverse_complement(self):
        pass 
    


