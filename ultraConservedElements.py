import multiprocessing
import numpy as np
import matplotlib.pyplot as plt
import os
import concurrent.futures
from multiprocessing import Process, Manager
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import random
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import RNA
import sys
from pkg_resources import cleanup_resources


from seqfold import fold



import matplotlib.pyplot as plt

# def run_blast_parallel(query):
#     return NCBIWWW.qblast(blast_program, database, query, expect=evalue, word_size=len(query_sequence))



def visualize_rna_structure(sequence, dot_bracket_notation, output_filename):

    # Create a new figure for the plot
    fig, ax = plt.subplots()

    # Plot the RNA secondary structure using Matplotlib
    ax.plot(range(1, len(sequence) + 1), [1 if c == '(' else -1 for c in dot_bracket_notation], color='blue', marker='o', linestyle='-', markersize=8)
    ax.set_ylim(-1.5, 1.5)
    ax.set_title('RNA Secondary Structure')

    # Save the plot to the desired output file
    plt.savefig(output_filename)










# def run_blast(query_sequence, blast_program="blastn", database="nt", evalue=1e-3, word_size=20, hitlist_size=10):
#     try:
#         # Perform BLAST search
#         result_handle = NCBIWWW.qblast(blast_program, database, query_sequence, expect=evalue, word_size=word_size, hitlist_size=hitlist_size, url_base='http://host.my.cloud.service.provider.com/cgi-bin/blast.cgi', format_object='Alignment')

#         print("Blast search successful.")

#         # Parse BLAST result
#         blast_records = NCBIXML.read(result_handle)
#         print("Dasfadsfads")
#         result_handle.close()

#         # Process and accumulate results in a string
#         result_string = ""
#         for alignment in blast_records.alignments:
#             for hsp in alignment.hsps:
#                 result_string += f"Alignment: {alignment.title}\n"
#                 result_string += f"Length: {alignment.length}\n"
#                 result_string += f"E-value: {hsp.expect}\n"
#                 result_string += f"{hsp.query}\n"
#                 result_string += f"{hsp.match}\n"
#                 result_string += f"{hsp.sbjct}\n\n"
#         print("adfsadsf")
#         return result_string

    # except Exception as e:
    #     # Handle exceptions here
    #     print(f"An error occurred during BLAST search: {str(e)}")
    #     return None

def findConservedElementsProcess(start, end, sequences, numBases, refSeq, refSeqIndex, conservedSequenceNumThreshold, numberOfKmersToFormFromRef, result_list):
    ultraConservedSequencesLocal = []  # Local list for each process
    posInRef = start

    while posInRef < end:  # goes through the kmers of the ref sequence
        ultraConservedSequence = []  # list of the indexes of sequences that were ultra conserved for a given numBases kmer
        window = refSeq[posInRef: posInRef + numBases]  # gets the ref sequence window to look for conserved elements

        for i in range(len(sequences)):
            if i != refSeqIndex:  # skips over the reference sequence
                curSequenceSection = sequences[i]  # current sequence

                if refSeq[posInRef:posInRef + numBases] in curSequenceSection:
                    if "N" not in refSeq[posInRef:posInRef + numBases]: #throws out any conserved element with missing values
                        ultraConservedSequence.append((i, curSequenceSection.index(refSeq[posInRef:posInRef + numBases]), posInRef))

        if len(ultraConservedSequence) >= conservedSequenceNumThreshold:
            ultraConservedSequencesLocal.append((ultraConservedSequence, window))

        posInRef = posInRef + 1  # updates the posInRef

    result_list.extend(ultraConservedSequencesLocal)

def findConserveredElements(sequences, numBases, refSeq, refSeqIndex, conservedSequenceNumThreshold, numberOfKmersToFormFromRef):
    ultraConservedSequences = []  # Store local results from each process

    num_processes = 282  # You can adjust the number of processes
    chunk_size = numberOfKmersToFormFromRef // num_processes

    with Manager() as manager:
        result_list = manager.list()  # Shared list to store results

        processes = []
        for i in range(num_processes):
            start = i * chunk_size
            end = start + chunk_size

            p = Process(target=findConservedElementsProcess,
                        args=(start, end, sequences, numBases, refSeq, refSeqIndex, conservedSequenceNumThreshold, numberOfKmersToFormFromRef, result_list))
            processes.append(p)
            p.start()

        for p in processes:
            p.join()

        # Concatenate the local results into the global list
        ultraConservedSequences.extend(result_list)

    return ultraConservedSequences

# Rest of your code remains unchanged


# Define a function for each process to handle the entire workflow
def process_data(parameters):
    # Set parameters for this process
    sequences, numBases, refSeq, refSeqIndex, conservedSequenceNumThreshold, numberOfKmersToFormFromRef, sequenceNames = parameters

    # Call findConservedElementsProcess
    conservedElements = findConserveredElements(sequences, numBases, refSeq, refSeqIndex, conservedSequenceNumThreshold, numberOfKmersToFormFromRef)
    #create top folder
    top_folder_name = f"conservedElements_refGenome{refSeqIndex}"
    os.makedirs(top_folder_name, exist_ok=True)
    # Create a new folder
    top_folder_name = f"conservedElements_refGenome{refSeqIndex}"
    subfolder_name =  f"conservedElements_{numberOfKmersToFormFromRef}mer_{numBases}bases_refGenome{refSeqIndex}"
    folder_name = os.path.join(top_folder_name, subfolder_name)
  
    os.makedirs(folder_name, exist_ok=True)

    # Write to file
    with open(os.path.join(folder_name, f"conservedElements.txt"), 'w') as out:
        out.write("Genomes: \n")
        out.write("------------------------------------------------------------------------\n")
        for i in range(len(sequenceNames)):
            out.write(f"sequence {i}: {sequenceNames[i]}\n")
        out.write("------------------------------------------------------------------------\n")
        out.write("Parameters: \n")
        out.write("----------------------------------------------------------\n")
        out.write(f"reference sequence: {refSeqIndex}\nnumber of sequences conserved threshold: {conservedSequenceNumThreshold}\n")
        out.write(f"number of bases per conserved element: {numBases}\nnumber of kmers to explore: {numberOfKmersToFormFromRef}\n")
        out.write("\n----------------------------------------------------------\n")
        out.write("Ultra Conserved Elements Found: \n")
        out.write("-----------------------------------------------------------------------------------------------------------\n")
        for i in range(len(conservedElements)):
            out.write(f"{conservedElements[i]}\n")
        out.write("-----------------------------------------------------------------------------------------------------------\n")

    # Create csv file with info to make a table
    with open(os.path.join(folder_name, f"conservedElementsTable.csv"), 'w') as out:
        out.write(f"Reference Genome: sequence {refSeqIndex}\n")
        for i in range(len(sequenceNames)):
            out.write(f"sequence {i}: {sequenceNames[i]}\n")
        
        out.write("Ultra Conserved Elements\n")
        for i in range(len(conservedElements)):
            element = conservedElements[i]
            sequenceConserved = element[1]
            sequencesWithElement = element[0]
            out.write(sequenceConserved + "\n,," + "Position in sequence,Position in Reference Sequence\n")
            
            for j in range(len(sequencesWithElement)):
                sequenceAtSeq = sequencesWithElement[j]
                out.write(f",sequence {sequenceAtSeq[0]},{sequenceAtSeq[1]},{sequenceAtSeq[2]}\n")

    

    x_labels = []
    percentages = []
    # Create bar graph
    plt.title("Conserved Elements")
    for i  in range(len(conservedElements)):
         element = conservedElements[i]
         sequenceInfo = element[1]
         sequenceNum= element[0]
         x_labels.append(str(sequenceInfo))
         percentages.append((len(sequenceNum)/ len(sequences) ) * 100)

    plt.ylabel("Percentage of Sequences Element is Conserved In (%)")

    plt.xlabel("Element Conserved")
   
 
    if len(conservedElements) < len(sequences):
        plt.figure(figsize=(len(x_labels) * 0.3,  7))
    else:
        plt.figure(figsize=(len(x_labels) * 0.3, 10))

    plt.bar(x_labels, percentages, width=.6)
    
    plt.ylabel("Percentage of Sequences Element is Conserved In (%)")

    plt.xlabel("Element Conserved")
     # Adjust x-axis label rotation and bottom margin
    plt.xticks(rotation=90, ha="right", fontsize=7)  # Rotate labels 45 degrees, set fontsize
    
    plt.subplots_adjust(bottom=0.3)  # Increase bottom margin
    plt.tight_layout()
  



    plt.savefig(os.path.join(folder_name, f"conservedElements.png"))
    plt.close()

    file_nameFNA = str(numberOfKmersToFormFromRef)+"mer_"+str(numBases)+"bases_refGenome_" + str(refSeqIndex) + "_conservedElementsExplicit_ForBlast.fna"

  # Write additional data to file
    with open(os.path.join(folder_name, file_nameFNA), 'w') as out:
        id = 0
        for label in x_labels:
            out.write(">Element"+ str(id) + "\n"  + label + "\n")
            id = id + 1
    
    newfolder= os.path.join(folder_name, "BLAST")
    os.makedirs(newfolder, exist_ok=True)

    # print(folder_name + str(x_labels[i]),f"_BLAST.txt")
    # for i in range(len(x_labels)):
        #  print(x_labels[i])

#         batch_size = len(x_labels)
#         for i in range(0, len(queries), batch_size):
#             batch_queries = queries[i]
#             for query in batch_queries:
#                 result_handle = NCBIWWW.qblast(
#                 blast_program="blastn",
#                 database="nt",  # Assuming you are using the nucleotide database
#                 query_sequence="YOUR_QUERY_SEQUENCE",  # Replace with your actual query sequence
#                 expect=1e-3,
#                 word_size=20,
#                 gapcosts="5 2",  # Adjust gap costs as needed
#                 filter=False,  # Disable low complexity filter
#                 megablast=True,  # Use megablast for highly similar sequences
#                 hitlist_size=10,  # Limit the number of hits returned
# )

    #     blast_result = run_blast(x_labels[i], word_size=len(x_labels[i]), hitlist_size=10)
    #     blast_result = blast_result + "fasfas"
    #     print(blast_result)
    #     #  print("adsfadsfds" + blastResult)
    #     #  print(x_labels[i])
    #     #  print(folder_name + str(x_labels[i]),f"_BLAST.txt")
    #     with open((folder_name + "/" + "BLAST/" + str(x_labels[i])+"_BLAST.txt"), 'w') as out:
    #          out.write(blast_result)

    # print(folder_name)
    # print(x_labels)
    #  #RNA fold file to see how well elements fold
  




          
    newfolder= os.path.join(folder_name, "RNA_folding_structure")
    os.makedirs(newfolder, exist_ok=True)
    for i in range(len(x_labels)):

        directory = folder_name + "/RNA_folding_structure/"+ str(x_labels[i]) +  "_RNA_Structure_Prediction.png"
    

            # Predict secondary structure using RNAfold
        rna = x_labels[i]
        rna = rna.replace('T', 'U')
        # Predict secondary structure using RNAfold
        

        structure, energy = RNA.fold(rna) #folds the rna returns dot structure
        



        visualize_rna_structure(rna,structure,directory) #saves the secondary rna stucture


  
   
    
    
  




    
def main():
  

    #opens our sequence files



    
   
    sequences = []
    sequenceNames = []



    folder_path = "/home/nwaka013/Desktop/CSCI 5481 Final Project - Ultra Conserved elements - Copy/Bacteria Genomes" # Replace this with the path to your folder

    # List all files in the folder
    files = os.listdir(folder_path)

    # Iterate through the files
    for file_name in files:
        # Construct the full path to the file
        file_path = os.path.join(folder_path, file_name)
        currSeq =  ""

        # Check if it's a file (not a subdirectory)
        if os.path.isfile(file_path):
            # Open and read the contents of the file
            with open(file_path, 'r') as file:
                first = True
                currSeq =  ""
                for line in file:
                     
                        if first == True:
                            first = False #makes it so that first < for first seq is ignored
                            name = line
                            name = name.replace(">", "")
                            name = name.replace("\n", "")
                            sequenceNames.append(name) #appends sequence name to sqeuence name
                        elif line[0] == ">":
                          
                            
                                currSeq =  currSeq + "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
                    
                        else:
                            
                            currSeq = currSeq + line #adds line to currseq
        sequences.append(currSeq.replace("\n", "")) #appends currSequence to list of sequences     
                          

    sequenceNames = [name.split(',')[0].strip() for name in  sequenceNames] #removes all the un needed information in the squence/genome name
    sequenceNames = [name.split('scaffold')[0].strip() for name in  sequenceNames] #removes all the un needed information in the squence/genome name


  
    
        

    file.close() #closes file after read



    # with open("s.txt", 'w') as out:
    #     for i in range(len(sequences)):
            
    #         out.write(sequences[i] + "\n") #writes genetic distance matrix to file


    # with open("n.txt", 'w') as out:
    #     for i in range(len(sequenceNames)):
            
    #         out.write(sequenceNames[i] + "\n") #writes genetic distance matrix to file

    #Tuning Parameters
    
    numBases = 15
    refSeq = sequences[19]   
    numSeqThreshhold = 1
    refSeqIndex = 19
    numberOfKmersToForm = 10000
    parameters_list = []
    num_processes = 40  # Adjust as needed

    # for i in range(): # setting parameters
    refSeq = sequences[6]
    refSeqIndex = 6
    for j in range(0, 41, 10):
        k = 1000
        numBases =  j
        # while k < 100000:
        #     numberOfKmersToForm = k
        parameters_list.append((sequences, numBases,refSeq, refSeqIndex, numSeqThreshhold, numberOfKmersToForm, sequenceNames)) #appends parameter tuple to list
            # k = k * 10
        max_workers = 1  # Set to 1 to ensure only one process at a time
        with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = [executor.submit(process_data, params) for params in parameters_list]

            # Wait for each task to complete before moving on to the next one
            for future in concurrent.futures.as_completed(futures):
                try:
                    result = future.result()
                    # Optionally process the result if needed
                except Exception as e:
                    print(f"An error occurred: {str(e)}")
            cleanup_resources()
        

    # conservedElements = findConserveredElements(sequences, numBases,refSeq, refSeqIndex, numSeqThreshhold, numberOfKmersToForm)





    # with open("conservedElements"+"_"+str(numberOfKmersToForm)+"mer_"+str(numBases)+"bases_"+"refGenome"+ str(refSeqIndex)+".txt", 'w') as out:
    #     out.write("Genomes: \n")
    #     out.write("------------------------------------------------------------------------\n")
    #     for i in range(len(sequenceNames)):
    #          out.write("sequence " + str(i) + ": " + str(sequenceNames[i]) +  "\n")
    #     out.write("------------------------------------------------------------------------\n")
    #     out.write("Parameters: \n")
    #     out.write("----------------------------------------------------------\n")
    #     out.write("reference sequence: "  + str(refSeqIndex) + "\n" +  "number of sequences conserved threshold: " +  str(numSeqThreshhold)  + "\n" +  "number of bases per conserved element: " + str(numBases) + "\n" + "number of kmers to explore: " + str(numberOfKmersToForm))#writes genetic distance matrix to file
    #     out.write("\n")
    #     out.write("----------------------------------------------------------\n")
    #     out.write("Ultra Conserved Elements Found: \n")
    #     out.write("-----------------------------------------------------------------------------------------------------------\n")
    #     for i in range(len(conservedElements)):
          
    #         out.write(str(conservedElements[i]) + "\n")
    #     out.write("-----------------------------------------------------------------------------------------------------------\n")


    # #create csv file with info to make table:
    # with open("conservedElementsTable"+"_"+str(numberOfKmersToForm)+"mer_"+str(numBases)+"bases_"+"refGenome"+ str(refSeqIndex)+".csv", 'w') as out:
    #      out.write("Reference Genome: " + "sequence " +  str(refSeqIndex) + "\n")
    #      for i in range(len(sequenceNames)):
    #           out.write("sequence " + str(i) + ": " + str(sequenceNames[i]) +  "\n")
         
    #      out.write("Ultra Conserved Elements\n")
    #      for i in range(len(conservedElements)):
    #           element = conservedElements[i]
    #           sequenceConserved = element[1]
              
    #           sequencesWithElement = element[0]
             
    #           out.write(sequenceConserved)
    #           out.write("\n,,"+ "Posititon in sequence" + "," + "Position in Reference Sequence" "\n")
              
    #           for j in range(len(sequencesWithElement)):
    #                 sequenceAtSeq = sequencesWithElement[j]
    #                 out.write(",sequence " + str(sequenceAtSeq[0]) + "," + str(sequenceAtSeq[1]) + "," + str(sequenceAtSeq [2]) + "\n")
     

         
   
    # plt.title("Conserved Elements")

    # x_labels = []
    # percentages = []
    # for i  in range(len(conservedElements)):
    #      element = conservedElements[i]
    #      sequenceInfo = element[1]
    #      sequenceNum= element[0]
    #      x_labels.append(str(sequenceInfo))
    #      percentages.append((len(sequenceNum)/ len(sequences) ) * 100)


    # if(len(conservedElements) < len(sequences)):
    #      plt.figure(figsize=(len(x_labels) * 0.3,  6))  # Adjust the width and height as needed
    # else:
    #     plt.figure(figsize=(len(x_labels) * 0.3,  len(sequences)))  # Adjust the width and height as needed
         
         
  
    # plt.ylabel("Percentage of Sequences Element is Conserved In (%)")

    # plt.xlabel("Element Conserved")
 

    # plt.bar(x_labels, percentages, width=.6)
    #  # Adjust x-axis label rotation and bottom margin
    # plt.xticks(rotation=45, ha="right", fontsize=7)  # Rotate labels 45 degrees, set fontsize
    
    # plt.subplots_adjust(bottom=0.3)  # Increase bottom margin
    # plt.xticks(rotation = 90)
    # plt.tight_layout()


    # plt.savefig("conservedElements"+"_"+str(numberOfKmersToForm)+"mer_"+str(numBases)+"bases_"+"refGenome"+ str(refSeqIndex)+".png")

    # with open("conservedElementsExplicit"+"_"+str(numberOfKmersToForm)+"mer_"+str(numBases)+"bases_"+"refGenome"+ str(refSeqIndex)+".txt", 'w') as out:
        
    #     for i in range(len(x_labels)):
    #         out.write(x_labels[i] +  "\n")

     
              










   


   



   

   
if __name__ == "__main__":
    
    main()
    
            
        




