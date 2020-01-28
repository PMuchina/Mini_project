#! /usr/bin/env python3
#
# python script for handling PDB files
# this program allows you to download PDB files(read-in) from the PDB database
# the read-in files are stored in a the user's working directory
# with the PDB files, the user can then choose from the options what he/she wants to extract from the file
# the user may opt to just view the contents or print the content in a file
#
# Peter Muchina 28/Jan/2020
# 

#this a dictionary to convert the three letter protein code to one letter
three2one = {
    'ALA' : 'A',    
    'ARG' : 'R',    
    'ASN' : 'N',    
    'ASP' : 'D',    
    'CYS' : 'C',    
    'GLU' : 'E',    
    'GLN' : 'Q',    
    'GLY' : 'G',    
    'HIS' : 'H',    
    'ILE' : 'I',    
    'LEU' : 'L',    
    'LYS' : 'K',    
    'MET' : 'M',    
    'PHE' : 'F',    
    'PRO' : 'P',    
    'SER' : 'S',    
    'THR' : 'T',    
    'TRP' : 'W',    
    'TYR' : 'Y',   
    'VAL' : 'V'
}

#this is an introduction interface that just aquints you to the program
name=input('Please enter your name\n:')
print(f'Welcome aboard {name}')
print(f'Hope you are ready for an awesome experience\n')
print('This program allows you to parse PDB files and extract as much information as you may require')

print('Before commencing, please enter a path to your working directory below')
print('eg../home/muchina/Python/python-mini-project-PMuchina/top20/')
print('Every file read-in will be stored here\n')

#this sets your working 
user_input = input("Enter the path to your working directory: ")


import sys
import os
import Bio

def read():   
    import sys
    import os
    import urllib.request
    
    assert os.path.exists(user_input)

    src = "http://www.rcsb.org/pdb/downloadFile.do" \
        + "?fileFormat=pdb&compression=NO&structureId=%s"

    #The pdb files you want to read are recorded in a text file(PDBfiles.txt)
    with open("PDBfiles.txt","w") as fh:
        filename= input("Enter PDB filename to be readin:\n")
        filename=filename.upper()
        fh.write(filename)

    #The pdb files are readin and stored in the PDBfiles directory    
    with open("PDBfiles.txt","r") as fh:
        for line in fh:
            pdb = line

        filename = os.path.join(user_input, pdb+".pdb")

        if os.path.isfile(filename) :
                print ("Already read in " + pdb)      
        else :
            print ("Downloading " + pdb)
            print("Please wait")
            urllib.request.urlretrieve( src % pdb, filename)
            print ("Successfully read in")
    
            


def search():
    import re
    Entries=os.listdir(user_input)
    print('Files already read in:')
    for entry in Entries:
        print(entry)

    filename =  input('Choose one file to work on:\n')
    raw_sequence=[]
    seff=[]  

    with open(user_input+filename,'r') as f:
        for line in f:
            if line.startswith("SEQRES"):
                list=(str.split(line[19:70]))
                for res in list:
                    seff.append(res)

                    raw_sequence.append(three2one[res])
    sequence=str(''.join(raw_sequence).replace(' ',''))
    pattern = re.compile(r'[GYL]A[PFW][TLVGMT]') #Find patterns starting with G or Y or L then A then P or F or W
    matches = pattern.finditer(sequence)
    for match in matches:
        print(match)
            
def write():
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Alphabet import IUPAC
    Entries=os.listdir(user_input)
    print('Files already read in:')
    for entry in Entries:
        print(entry)
    filename =  input('Choose one file to work on:\n')
    
    options = input("""
    1-Write out coordinate file
    2-Write out sequence (Fasta format)
    
    """)
    
    if options== "1":
        sub_options=input("""
        a- All atom
        b- Backbone atoms
        c-Alpha carbon atoms only
        """)
        
        if sub_options == "a":
            input_file= input("Please enter a .txt file to write the coordinates:\n")
            with open(user_input+filename,'r') as f:
                with open(input_file,"w") as fh:
                    for line in f:
                        if line.startswith("ATOM"):
                            fh.write(line)
                    print(f'Coordinates written succesfully in: {input_file}')

        if sub_options == "b":
            input_file= input("Please enter a .txt file to write the coordinates:\n")
            with open(user_input+filename,'r') as f:
                with open(input_file,"w") as fh:
                    for line in f:
                        if line[:4] == 'ATOM' and line[12:16] == " CA ":
                            fh.write(line)
                        elif line[:4] == 'ATOM' and line[12:15] == " N ":
                            fh.write(line)
                        elif line[:4] == 'ATOM' and line[12:15] == " C ":
                            fh.write(line)
                        elif line[:4] == 'ATOM' and line[12:15] == " O ":
                            fh.write(line)
                    print(f'Coordinates written succesfully in: {input_file}')

        if sub_options == "c":
            input_file= input("Please enter a .txt file to write the coordinates:\n")
            with open(user_input+filename,'r') as f:
                with open(input_file,"w") as fh:
                    for line in f:
                        if line[:4] == 'ATOM' and line[12:16] == " CA ":
                            fh.write(line)
                    print(f'Coordinates written succesfully in: {input_file}')
                        
    if options== "2":
         sub_options=input("""
    1-SEQRES Sequence
    2-Coordinate Sequence
    3-Alignment Sequence
    """)
    
    
    if sub_options == "1":
        input_file= input("Please enter a .fasta file to write the SEQRES sequence:\n")
        input_file= input("Please enter a .fasta file to write the SEQRES sequence:\n")
        with open(user_input+filename,'r') as f:
            with open(input_file,"w") as fh:
                print("Converting to FASTA...")
                print("Just smile...")
                print("Please wait...")
                for strLine in f:
                    if strLine.startswith("SEQRES"):
                        line=strLine[17:]
                        strLine=line.rstrip("\n")
                        record = SeqRecord(Seq(strLine,IUPAC.protein),id="PK_025292.1", name=filename, 
                                           description=f'PDB files in {user_input}')
                        fh.write(record.format('fasta'))

                print("Written succesfully")
                    
    if sub_options == "2":
        input_file= input("Please enter a .fasta file to write the Coordinate sequence:\n")    
        with open(user_input+filename,'r') as f:
            with open(input_file,"w") as fh:
                print("Converting to FASTA...")
                print("Just smile...")
                print("Please wait...")
                for strLine in f:
                    if strLine.startswith("ATOM"):
                        line=strLine[30:-27]
                        strLine=line.rstrip("\n")
                        record = SeqRecord(Seq(strLine,IUPAC.protein),id="PK_025292.1", name=filename, 
                                           description=f'PDB files in {user_input}')
                        fh.write(record.format('fasta'))

                print("Written succesfully")
                
    if sub_options == "3":
        input_file= input("Please enter a .txt file to write the Alignment sequence:\n") 
        raw_sequence=[]
        seff=[]  

        with open(user_input+filename,'r') as f:
            with open(input_file,"w") as fh:
                for line in f:
                    if line.startswith("SEQRES"):
                        list=(str.split(line[19:70]))
                        for res in list:
                            seff.append(res)

                            raw_sequence.append(three2one[res])
                sequence=str(''.join(raw_sequence).replace(' ',''))
                record = SeqRecord(Seq(sequence,IUPAC.protein),id="PK_025292.1", name=filename, 
                                           description=f'PDB files in {user_input}')
                fh.write(record.format('fasta'))

                print("Written succesfully")

def information() :   
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Alphabet import IUPAC
    Entries=os.listdir(user_input)
    print('Files already read in:')
    for entry in Entries:
        print(entry)
    filename =  input('Choose one file to work on:\n')

    sub_options=input("""
    1-Display coordinate sequence (as explained above)
    2-Display SEQRES sequence (as explained above)
    3-Display Alignment sequence (as explained above)
    4-Display all non-water ligands in the protein (if any)

    """)

    if sub_options=="1":
        with open(user_input+filename,'r') as f:
            for line in f:
                if line.startswith("ATOM"):
                    print(line[30:-27])

    if sub_options=="2":
        with open(user_input+filename,'r') as f:
            for line in f:
                if line.startswith("SEQRES"):
                    line=line[17:].rstrip()
                    print(line)

    if sub_options=="3":
        raw_sequence=[]
        seff=[]  
        with open(user_input+filename,'r') as f:
            for line in f:
                if line.startswith("SEQRES"):
                    list=(str.split(line[19:70]))
                    for res in list:
                        seff.append(res)

                        raw_sequence.append(three2one[res])
            sequence=str(''.join(raw_sequence).replace(' ',''))
            record = SeqRecord(Seq(sequence,IUPAC.protein),id="PK_025292.1", name=filename, 
                                       description=f'PDB files in {user_input}')
            print(record)

def alignment():
    raw_sequence1=[]
    seff1=[] 
    raw_sequence2=[]
    seff2=[]
    
    Entries=os.listdir(user_input)
    print('Files already read in:')
    for entry in Entries:
        print(entry)
    print("Please select two PDB files to be aligned")
    
    filename1=input("Enter the name of the first file:\n")
    filename2=input("Enter the name of the second file:\n")
    output_file=input("Enter a .txt file to write the alignment to:\n")
    
    with open(user_input+filename1,'r') as f:
        for line in f:
            if line.startswith("SEQRES"):
                list=(str.split(line[19:70]))
                for res in list:
                    seff1.append(res)

                    raw_sequence1.append(three2one[res])
    sequence1=str(''.join(raw_sequence1).replace(' ',''))
    #print(sequence1)
    X=sequence1[:50]
                   
                
    with open(user_input+filename2,'r') as f:
        for line in f:
            if line.startswith("SEQRES"):
                list=(str.split(line[19:70]))
                for res in list:
                    seff2.append(res)

                    raw_sequence2.append(three2one[res])
    sequence2=str(''.join(raw_sequence2).replace(' ',''))
    #print(sequence2)
    Y=sequence2[:50]
                               
    # Import pairwise2 module
    from Bio import pairwise2

    # Import format_alignment method
    from Bio.pairwise2 import format_alignment

    # Define two sequences to be aligned

    #X 
    #Y

    # Get a list of the global alignments between the two sequences ACGGGT and ACG satisfying the given scoring
    # A match score is the score of identical chars, else mismatch score.
    # Same open and extend gap penalties for both sequences.
    alignments = pairwise2.align.globalms(X, Y, 2, -1, -0.5, -0.1)

    # Use format_alignment method to format the alignments in the list
    with open(output_file,"w") as fh:
        for a in alignments:
            fh.write(format_alignment(*a))
    print(f'Alignment Written succesfully to {output_file}')


def help():
    print("These are the various options in the programme:")
    print("""
    R – Read in PDB files
    S – Search option
    W – Write out option
    I – Information option
    A – Alignment option
    H – Help option
    Q – quit

    """)
    
    
    print("""
    1.The read option allows you do download any PDB file which is stored in the system to be used later on
    
    2.The search option allow you to find all the gycosylation sites in the PDB file
    
    3.The write option allow you to parse through the file and write the desired sequences into files of
    your choice
    
    4.The information option displays the desired sequences without writing them to files
    
    5.The alignment option allows you to extract sequences from two PDB files, align them and write the 
    results into a file of your choice
    
    6.The help options displays for you all the options for easier navigation through the program
    
    7.The quit option terminates the program
    
    
    """)
    
    
def quit():
    print("Goodbye....")
    os._exit


def logo():
    print(f'\033[1;32;40m{"-"* 100}')
    print(f'\033[1;32;40m{"-"* 100}\n')
    result_str="";    
    for row in range(0,7):    
        for column in range(0,7):     
            if (column == 1 or ((row == 0 or row == 3) and column > 0 and column < 5) or ((column == 5 or column == 1) and (row == 1 or row == 2))):  
                result_str=result_str+"p"    
            else:      
                result_str=result_str+" "    
        result_str=result_str+"\n"    
    print(f'\033[1;32;40m{result_str}');

    result_str="";    
    for row in range(0,7):    
        for column in range(0,7):     
            if (column == 1 or ((row == 0 or row == 6) and (column > 1 and column < 5)) or (column == 5 and row != 0 and row != 6)):  
                result_str=result_str+"D"    
            else:      
                result_str=result_str+" "    
        result_str=result_str+"\n"    
    print(f'\033[1;32;40m{result_str}');

    result_str="";    
    for row in range(0,7):    
        for column in range(0,7):     
            if (column == 1 or ((row == 0 or row == 6 or row ==3) and (column > 1 and column < 5)) or (column == 5 and row != 0 and row != 6)):  
                result_str=result_str+"B"    
            else:      
                result_str=result_str+" "    
        result_str=result_str+"\n"    
    print(f'\033[1;32;40m{result_str}');

    print("-"* 100)
    print("-"* 100)   
    
        


    
    def menu():
        option=input("""
        R – Read in PDB files
        S – Search option
        W – Write out option
        I – Information option
        A – Alignment option
        H – Help option
        Q – quit

        """)
        if option !='Q':
            if option.upper() in ('R','S','W','I','A','H'):
                if option.upper() == 'R':
                    read()
                    menu()
                if option.upper() == 'S':
                    search()
                    menu()
                if option.upper() == 'W':
                    write()
                    menu()
                if option.upper() == 'I':
                    information()
                    menu()
                if option.upper() == 'A':
                    alignment()
                    menu()
                if option.upper() == 'H':
                    help()
                    menu()
                if option.upper() == 'Q':
                    quit()
    menu()
logo()



