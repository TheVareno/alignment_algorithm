#/opt/anaconda3/bin/python

import tkinter as tk 
from tkinter import *
from alignment.sequence import Sequence
from alignment.alignment import Alignment 
from alignment.algorithms import Algorithms
import re

def align_sequences():
    # get the data from GUI
    # input 1 insereted by user
    sequence1 = input_seq1.get() 
    # input 2 inserted by user
    sequence2 = input_seq2.get()
    
    # aligment type selected by user
    align_type = ''
    
    # create sequence pattern to identify the sequence type 
    nt_pattern = re.compile(r"^[ACTG]+$")
    # check the sequence type
    seq_type = ''
    if nt_pattern.match(sequence1):
        seq_type = 'nt'
    else:
        seq_type = 'aa' 
    
    # create sequence object     
    seq1 = Sequence(sequence1, seq_type)
    seq2 = Sequence(sequence2, seq_type)    
    
    # run the algorithm using Algorithm object 
    algorithm = Algorithms(align_type, seq_type)
    if align_type == 'global': 
        alignment = algorithm.global_alignment(seq1, seq2)
    elif align_type == 'local':
        alignment = algorithm.local_alignment(seq1, seq2)
    else : 
        pass 
    
    # the first output is an alignment object produced by algorithm object 
    return alignment
        

def update_result_label():
    # get the alignment object from align_sequence() func
    alignment = align_sequences()   
    # Update the result label with alignment objects properties 
    
    # add the third widgets : labels to show the result 
    result_label_scorecost = tk.Label(main_window, text=alignment.alignment_score)
    result_label_scorecost.pack()
    result_label_alignment = tk.Label(main_window, text=str(alignment.aligned_sequences))
    result_label_alignment.pack()
    

if __name__ == "__main__" : 
    # create ui fundamental
    main_window = tk.Tk()
    # adjust the main window size 
    main_window.geometry('500x500')

    # add the first widget : two input fields to get the sequence data 
    input_seq1 = tk.Entry(main_window)  
    input_seq1.pack()
    input_seq2 = tk.Entry(main_window)
    input_seq2.pack()


    # add the second widget : a button to initilize the alignment
    align_btn = tk.Button(main_window, text='Start the Alignmentt', command=update_result_label)
    align_btn.pack() 

    # add the dropdown menu to choose alignment type
    alignment_type_options = ["global", "local"]
    clicked = StringVar()
    clicked.set("global")
    drop = OptionMenu(main_window, clicked, *alignment_type_options)
    drop.pack()

    main_window.mainloop()