#!/home/genouest/dyliss/norobert/miniconda3/envs/prolific/bin python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 10 2023

@author: ytirlet, norobert
"""

from __future__ import print_function
import os,sys
import os.path
from optparse import OptionParser
import pandas as pd
import numpy as np
import shutil
import gzip
import re


# FUNCTIONS ---------------------------------------------------------------------------------

def rename(name) :
    """ from a [str] name, returns the equivalent str with genera shortened when listed, 
        and the annotation level converted to a 1-or-2-letter id """
    dico_prefix = {'Anaerobutyricum':'Ab','Anaerostipes':'As','Bifidobacterium':'B','Coprococcus':'Co','Clostridium':'C','Enterococcus':'E', 'Faecalibacterium':'Fa','Fructilactobacillus':'F','Lactobacillus':'Lb', 'Lacticaseibacillus':'Lb','Lactiplantibacillus':'Lb' ,'Limosilactobacillus':'Lb','Levilactobacillus':'Lb','Lentilactobacillus':'Lb','Latilactobacillus':'Lb','Schleiferilactobacillus':'Lb','L':'Lco', 'Lactococcus':'Lco',  'Leuconostoc':'Leu', 'Pediococcus':'Pe', 'Propionibacterium':'Pr', 'P':'Pr', 'Streptococcus':'St','Periweissella':'W','Weissella':'W'}
    dico_suffix = {'Complete':'C', 'complet':'C', 'C':'C', 'contig':'co','Contig':'co','scaffold':'S','Scaffold':'S', 'S':'S', 'Plasmid':'P', 'plasmide':'P', 'P':'P'}

    parts = name.split('_')
    prefix = parts[0]
    middle = parts[1:-1]
    suffix = parts[-1]
    # prefix
    if prefix in dico_prefix :
        new_name = dico_prefix[prefix] + '_'
    else :
        new_name = prefix + '_'
    for i in range(len(middle)) :
        new_name += middle[i] + '_'
    # suffix
    if suffix in dico_suffix :
        new_name += dico_suffix[suffix]
    else :
        new_name += suffix + '_S'

    return new_name

def forbidden(name) :
    """ in a given string, replace all forbidden characters ('.', ':', '__', '-_') by a '-' """
    while '.' in name or ':' in name or '__' in name or '-_' in name :
        name = name.replace('.','-')
        name = name.replace(':','-')
        name = name.replace('__','-')
        name = name.replace('-_','-')
    return name

def bigprint(message): 
    delimitation = "-------------------------------------------"
    print("\n{}\n{}\n{}\n".format(delimitation,message,delimitation))
    return

def move(source, dest) : 
    try :
        shutil.move(source,dest)
    except PermissionError:
        print("Some rights are missing to move {} to {}".format(source,dest))

def mkdir(path) : 
    if not os.path.exists(path):
        try :
            os.mkdir(path)
        except PermissionError:
            print("Some rights are missing to create {}".format(path))
        except Exception as e:
            print(f"An error occurred (mkdir): {e}")

def remove(list_path) : 
    for path in list_path : 
        if os.path.exists(path):
            os.remove(path)
        elif os.path.isdir(path):
            shutil.rmtree(path)

def decompress_gzip_file(file_path, name, suppr_zip):
    if file_path.endswith(".gz"):
        dir_path = file_path.split("/")
        dir_path = "/".join(dir for dir in dir_path[:-1])
        with gzip.open(file_path, 'rb') as f_in, open(dir_path+"/"+ name +".fasta", 'wb') as f_out:
            f_out.write(f_in.read())
        if suppr_zip == True : 
            os.remove(file_path)  

def decompression(new_names, path_to_all_data) :
        print("decompressing input files... ")
        for new_name in new_names :
            genome_dir = path_to_all_data+new_name

            ## unzipping the files if tared 
            for root, dirs, fastas in os.walk(genome_dir):
                for fasta in fastas :
                    if fasta.endswith(".gz"):
                        if re.search(new_name, genome_dir+fasta): ## if name contains directory's
                            ## False : don't delete archive
                            decompress_gzip_file(f"{genome_dir}/{fasta}", new_name, False)  
        print("decompressing over.")       

def missing_or_empty(file_path):
    return not os.path.exists(file_path) or os.stat(file_path).st_size == 0

def checkList(annots, annotools, path):
    """
    Input :
        annots : list of x (preferentially 3) sublists of found files 
        annotools : list of annotation tool's name to get index
        path : path of checked dir
    Output : 
        True if len(annot)!=3 or if no file is missing from any subdir 
        (i.e. if there's no fail in annotation)
        False if any file is missing 
    """
    if len(annots) == 3 : 
        for i, annot in enumerate(annots):
            possibilities = [0, 1, 2]
            possibilities.remove(i)
            for name in annot :
                for other in possibilities :
                    if name not in annots[other] :
                        message = f"ERROR : {name} from {annotools[i]} not found in {path}{annotools[other]}"
                        return False, message
        return True
    else : 
        print(f"{len(annots)} different annotation(s) instead of 3 expected. Verification function won't be launched.")
        return True

# ---------------------------------------------------------------------------------------------


def main() :
    parser = OptionParser()
    parser.add_option("-i", "--input", dest="input",help="Path to the folder where the genomes are")
    parser.add_option("-o", "--output", dest="output",help="Path to the folder where you want to put the results in")
    parser.add_option("--tax", dest="all_taxon",help="path of the all_taxon.tsv file")
    parser.add_option("--padmet_ref", dest="path_to_padmet_ref", help="Path to the reference database in Padmet format.")
    parser.add_option("--ptsc",dest="ptsc", help="Path to root folder.")
    parser.add_option("--ptsi",dest="ptsi", help="Name of the singularity image.")
    parser.add_option("--pwy",dest="pwy_fold", help="Path to the folder with the pathways.txt files for all wanted metabolites.")
    parser.add_option("--strain",dest="strain", help="Path to the strains file.")
    parser.add_option("--annot",dest="annot",help="Annotation tool(s) to use between 'prokka' (default), 'eggnog' and 'bakta'. If several annotation tools to use, write them comma-separated.")
    parser.add_option("--egg_path",dest="egg_path",help="Path to the eggnog database, mandatory if you want to use eggnog as annotation tool.")
    parser.add_option("--bak_path",dest="bak_path",help="Path to the bakta database, mandatory if you want to use bakta as annotation tool.")
    parser.add_option("-r","--rename",action="store_true",dest="rename", help="Renames of the strains with abreviations.")
    parser.add_option("-a","--asko", action="store_true", dest="asko", help="Launch the creation of the askomics files.")
    parser.add_option("-v","--verbose",action="store_true",dest="verbose", help="Activate verbose.")
    parser.add_option("-k","--keep_faa", action="store_true", dest="keep_faa", default=False, help="Keep .faa files that can be need to use other annotation software like eggNOG-mapper")
    parser.add_option("-u","--unzip", action="store_true", dest="unzip", default=False, help="write this flag if your files are gziped")
    parser.add_option("-c","--cpus", dest="cpus", default=20, help="Give the number of available CPUs")
    parser.add_option("-q", "--quick", action="store_true", dest="quick", help="Bypass most of the computation if results files are already generated")
    
    
    (options,args) = parser.parse_args()

    if len(args) >=1:
        print(parser.print_help())
        parser.error("Incorrect number of arguments")
        return
    
    path_to_all_data = options.input
    files = os.listdir(options.input)
    path_to_scratch = options.ptsc
    path_to_singularity = options.ptsi
    output_path = options.output
    all_taxon = options.all_taxon
    strain_file = options.strain
    nbThreads=options.cpus
    quick = options.quick

    if options.annot :
        annotation = options.annot.split(",")   # annotation becomes a list if different
    else :
        annotation = 'prokka'

    
    if not quick : 
        ## Creating output directory
        os.makedirs(output_path, exist_ok=True)

        #-------------------------------------------------------
            # UNCOMPRESS AND RENAME THE FILES
        #-------------------------------------------------------    
        ## unzipping the files if tared 
        if options.unzip :
            print("decompressing input files... ")
            for name in files : 
                genome_dir = path_to_all_data+name
                for root, dirs, fastas in os.walk(genome_dir):
                    for fasta in fastas :
                        if fasta.endswith(".gz"):
                            if re.search(name, genome_dir+fasta):
                                ## False : don't delete archive
                                decompress_gzip_file(f"{genome_dir}/{fasta}", False)  
            print("decompressing over.")

        new_names = []
        for name in files :
            if options.rename :
                new_name = forbidden(rename(name))
            else :
                new_name = forbidden(name)
            new_names.append(new_name)
            genome_dir = path_to_all_data+name

        if options.rename :
            print("Conversion of fna files to fasta (if any) and renaming input...")
        
            ## renaming files with their directory's name, without forbidden caracters
            for file in os.listdir(genome_dir):
                if file.endswith(".fna") or file.endswith(".fasta"):
                    move(genome_dir+"/"+file, genome_dir+"/"+new_name+".fasta")
                if file.endswith(".fasta"):
                    move(genome_dir+"/"+file, genome_dir+"/"+new_name+".fasta")
            
            ## renaming the directory
            if os.path.exists(genome_dir):
                move(genome_dir, path_to_all_data+new_name)
        print("done.")       

        #-------------------------------------------------------
            # ANNOTATION
        #-------------------------------------------------------

        if 'prokka' in annotation  :
            #-------------------------------------------------------
                # USING PROKKA FOR ANNOTATION
            #-------------------------------------------------------
            print("Prokka annotation launched.")
            mkdir(output_path + 'prokka')

            for new_name in new_names : 
                command_pro = f"prokka {path_to_all_data}{new_name}/* --outdir {output_path}prokka/{new_name} --prefix {new_name} --compliant --force --cpus {nbThreads}"
                bigprint(command_pro)
                os.system(command_pro)
                ## --compliant       Force Genbank/ENA/DDJB compliance
    
                prok_file = f"{output_path}prokka/{new_name}/{new_name}"
                if os.path.exists(prok_file+".gbf"):
                    move(prok_file+".gbf",prok_file+".gbk")     #transform .gbf to .gbk
                remove([f"{prok_file}.ecn",f"{prok_file}.err",f"{prok_file}.ffn",f"{prok_file}.fixed*",f"{prok_file}.fsa",f"{prok_file}.gff",f"{prok_file}.log",f"{prok_file}.sqn",f"{prok_file}.tbl",f"{prok_file}.val"])
                if options.keep_faa == False :
                    remove([prok_file+".faa "])
            
        if 'eggnog' in annotation :
            print("Eggnog annotation launched.")
            path_to_egg = options.egg_path
            #-------------------------------------------------------
                # USING EGGNOG-MAPPER FOR ANNOTATION
            #-------------------------------------------------------
            mkdir(output_path + 'eggnog')

            for new_name in new_names :
                mkdir(output_path + 'eggnog/' + new_name)
                command_egg = f"emapper.py -i {path_to_all_data}{new_name}/{new_name}.fasta -o {new_name} --cpu {nbThreads} --itype genome --data_dir {path_to_egg} --output_dir {output_path}eggnog/{new_name}/ --dbmem --genepred prodigal --override"
                bigprint(command_egg)
                os.system(command_egg)
                
                ## conversion of eggnog output to gbk
                genom = path_to_all_data + new_name + '/' + new_name + '.fasta'
                prot = output_path + 'eggnog/' + new_name + '/' + new_name + '.emapper.genepred.fasta'
                gff = output_path + 'eggnog/' + new_name + '/' + new_name + '.emapper.genepred.gff'
                annot = output_path + 'eggnog/' + new_name + '/' + new_name + '.emapper.annotations'
                out_file = output_path + 'eggnog/' + new_name + '/' + new_name + '.gbk'
                command_egg2gbk = f'emapper2gbk genomes -fn {genom} -fp {prot} -g {gff} -a {annot} -o {out_file} -gt eggnog -c {nbThreads}'
                # bigprint(command_egg2gbk)
                # os.system(command_egg2gbk)

                
        if 'bakta' in annotation :
            print("Bakta annotation launched.")
            path_to_bak = options.bak_path

            mkdir(output_path + 'bakta')

            for new_name in new_names :
                mkdir(output_path + 'bakta/' + new_name)
            
                command = f"bakta --db {path_to_bak} {path_to_all_data}{new_name}/{new_name}.fasta --output {output_path}/bakta/{new_name} --prefix {new_name} --compliant --force --threads {nbThreads}"
                ## --compliant      Force Genbank/ENA/DDJB compliance
                ## --force          Force overwriting existing output folder
                bigprint(command)
                os.system(command)

                ## removing unused files
                unused_files=[".embl", ".faa", ".ffn", ".fna", ".gff3", ".hypotheticals.faa", ".hypotheticals.ftsv", ".json", ".log", ".png", ".svg", ".tsv"]     
                for extension in unused_files : 
                    remove([f"{output_path}/bakta/{new_name}/{new_name}{extension}"])


        else :
            raise ValueError("The specified annotation tool is not recognized. Please retry with 'eggnog', 'bakta' or 'prokka'. Default is 'prokka'.")


        # -------------------------------------------------------
            # CREATE TAXON ID FILE
        # -------------------------------------------------------
        df_to_write = pd.DataFrame(columns = ["species", "taxon_id", "corresponding_file"])
        df_taxons = pd.read_csv(all_taxon, sep='\t')
        for index, row in df_taxons.iterrows() : 
            if options.rename :
                new_name = forbidden(rename(row["corresponding_file"]))    
            else : 
                new_name = forbidden(row["corresponding_file"])
            if new_name in new_names :
                df_to_write.loc[len(df_to_write)] = [new_name, row["taxon_id"], new_name]

        for annotool in annotation : 
            tax_file = output_path + annotool + '/taxon_id.tsv'
            df_to_write.to_csv(tax_file, sep="\t", index=False)  

        #-------------------------------------------------------
            # RUNNING MPWT USING SINGULARITY TO CREATE .dat FILES 
        #-------------------------------------------------------
        mkdir(output_path + 'mpwt')
        for annotool in annotation :
            mkdir(output_path + 'mpwt/' + annotool)

            ## checking if mpwt has successfully run before
            path =  os.path.join(output_path, "mpwt", annotool)
            dat_dirs = [d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d))]  ## lists subdirectories names
            dat_dirs = [d for d in dat_dirs if d.startswith("GCF")]                           ## filters for those which start with GCF
            print(f"{len(dat_dirs)} mpwt repositories found out of {len(new_names)} genomes to process")
            if len(dat_dirs) != len(new_names):
                command_mpwt = f"singularity exec -B {path_to_scratch}:{path_to_scratch} {path_to_scratch}{path_to_singularity} mpwt -f {output_path}{annotool}/ -o {output_path}mpwt/{annotool}/ --cpu {nbThreads} --patho --flat --clean --md -v"
                ## --patho : Launch PathoLogic inference on input folder
                ## --flat : Create BioPAX/attribute-value flat files
                ## --clean : Delete all PGDBs in ptools-local folder or only PGDB from input folder
                ## --md : Move the dat files into the output folder
                bigprint(command_mpwt)
                os.system(command_mpwt)
            else :
                print("Nothing to process here, moving on.")
    
        #-------------------------------------------------------
            # CONVERT .dat INTO .padmet FILES
        #-------------------------------------------------------
        path_to_padmet_ref= options.path_to_padmet_ref
        padmet_output = output_path + 'padmet'
        mkdir(padmet_output)

        ## checking if mpwt ran correctly
        mpwt_path = f"{output_path}mpwt/"
        list_dat_files = [os.listdir(f"{mpwt_path}{annotool}") for annotool in annotation]
        check = checkList(list_dat_files, annotation, mpwt_path)
        if type(check) is tuple: ## if number of files from subdirs of mpwt are not matching :
            bigprint(check[1])
            sys.exit(1) 
        
        for annotool in ["prokka"] : # annotation :
            dat_files = os.listdir(f"{output_path}mpwt/{annotool}")
            print(f"{annotool} : {len(dat_files)}")
            for name in new_names :
                if missing_or_empty(os.path.join(padmet_output, name, f"{name}_{annotool}.padmet")):
                    ## create files in commune directories for annotations of the same genome 
                    mkdir(f"{padmet_output}/{name}")
                    command_pgdb2padmet_source = f"singularity run -B {path_to_scratch}:{path_to_scratch} {path_to_scratch}{path_to_singularity} padmet pgdb_to_padmet --source=annot_{annotool} --pgdb={output_path}mpwt/{annotool}/{name}/ --output={padmet_output}/{name}/{name}_{annotool}.padmet --extract-gene --no-orphan --padmetRef={path_to_padmet_ref} -v"
                    bigprint(command_pgdb2padmet_source)
                    os.system(command_pgdb2padmet_source)

        # -------------------------------------------------------
            # MERGE .padmet FILES
        # -------------------------------------------------------            
        output_merged=output_path + 'merged_padmet/'
        mkdir(output_merged)
        for name in new_names :  
            if missing_or_empty(os.path.join(output_merged, name, ".padmet")):
                ## Check if 3 files are present for merging
                nb_of_padmets=len(os.listdir(os.path.join(padmet_output, name)))
                if nb_of_padmets == len(annotation) :
                    ## Merge annotation files for each genomes into one
                    command_padmet2padmet = f"singularity run {path_to_scratch}{path_to_singularity} padmet padmet_to_padmet --to_add={padmet_output}/{name}/ --output={output_merged}{name}.padmet -v"
                    bigprint(command_padmet2padmet)
                    os.system(command_padmet2padmet)
                else :
                    print(f"ERROR : Something wrong in {padmet_output}/{name}, couldn't merge padmets")
                    sys.exit()
    
        #-------------------------------------------------------
            # COMPARE THE .padmet FILES to get the tsv_files
        #------------------------------------------------------- 
        output_tsv = output_path + 'tsv_files'
        mkdir(output_tsv)
        command_compare_padmet = f"singularity run {path_to_scratch}{path_to_singularity} padmet compare_padmet --padmet={output_merged} --output={output_tsv} -v"
        bigprint(command_compare_padmet)
        os.system(command_compare_padmet)

    #-------------------------------------------------------
        # ANALYSE OF THE METABOLIC PATHWAYS
    #-------------------------------------------------------
    output_metabo = output_path + 'metabo_files'
    mkdir(output_metabo)
    reactions = output_path + 'tsv_files/reactions.tsv'
    p_files = os.listdir(options.pwy_fold)
    path_dir = options.pwy_fold
    nb_col = len(files) + 1
    all_not_found = []      ## exhaustive list of never-found reactions

    ## linking file names to genome names   
    strainfile_content = open(strain_file)
    filename2strain = {}
    filename2status = {}
    for line in strainfile_content :
        line = line.split('\t') 
        filename2strain[forbidden(line[2][:-1])] = line[0]
        if line[1].strip() == '' :                          ## trying to replace '' by "null"
            filename2status[forbidden(line[-1][:-1])] = "null"
        else :
            filename2status[forbidden(line[-1][:-1])] = line[1]

    ## for each metabolite, extracts the list of the pathway's reactions 
    for metabo in sorted(p_files) :
        output_file = output_path + 'metabo_files/' + metabo[:-4] + '.tsv'
        df = pd.read_csv(reactions, sep='\t', header=0, low_memory=False)
        
        ## get reactions list for the current metabolite 
        list_rxn = []
        fp = open(path_dir + metabo)
        for line in fp :
            if line[:-1] != "" : 
                line = line.replace(" ", "").replace("\n", "")
                list_rxn.append(line)

        ## extract reactions from reactions.tsv matching pathway's 
        df = df[df['reaction'].isin(list_rxn)]
        df = df.iloc[:,:nb_col]

        ## writing the tab in a .csv file with additionnal information
        df_t = df.T                     # transpose df
        rownames = df_t.index.values    # extract index columns
        rownames = list(rownames)
        tab = df_t.values               # df content into arrays
        column_names = tab[0]
        rows = np.array([rownames]).T   # column 1 : lines names

        ## shift names
        rows_strain = [['reaction']]
        for row in rows : 
            if "reaction" not in row :  
                rows_strain.append([row[0]])
        rows_strain = np.array(rows_strain)
        # tab = np.append(rows,tab,axis = 1)
        tab = np.append(rows_strain,tab,axis = 1)  
        sort_tab = tab[1:]
        sort_tab = sort_tab[sort_tab[:,0].argsort()]
        row_names = [tab[0,0]] + list(sort_tab[:,0])
        tab = sort_tab[:,1:]            # reactions presence/absence
        reaction_nb = len(list_rxn)

        # writing the tab in a csv file
        fo = open(output_file,"w")
        fo.write(row_names[0] + '\t')
        for name in column_names :
            fo.write(str(name) + '\t')
        
        ## get the reactions not found
        not_found = [react for react in list_rxn if react not in column_names]
        for react in list_rxn :
            if react not in column_names : 
                all_not_found.append(react)
        real_reaction_nb = len(list_rxn) - len(not_found)

        ## add them to the csv file
        for react in not_found :
            fo.write(react + '\t')
        fo.write('Number of possessed reactions\tTotal number of Reactions\tCompletion percent\tAdj Compl Pct\n')

        maxi = 0
        for i in range(1,len(row_names)) :
            react_count = 0
            fo.write(row_names[i] + '\t')
            for react_presence in tab[i-1] :
                react_count += int(react_presence)
                fo.write(str(react_presence) + '\t')
            # filling with zeros the columns of not-found-reactions
            for react in range(len(not_found)) :
                fo.write('0\t')
            # adding percentage of completion
            #print(f"{metabo} : {react_count} / {reaction_nb}")
            percent = round ((react_count / reaction_nb) * 100, 2)
            if real_reaction_nb == 0 : 
                adj_percent = 0
            else :
                adj_percent = round ((react_count / real_reaction_nb) * 100, 2)
            if react_count > maxi : 
                maxi = react_count
            fo.write(str(react_count) + '\t' + str(reaction_nb) + '\t' + str(percent) + "\t" + str(adj_percent) + '\n')
        fo.close()

        ## few informations about reactions repartition
        print(f"{metabo} : max. {maxi}/{reaction_nb} ; {len(not_found)} never found ({', '.join(reaction for reaction in not_found)})")
    
    print(f"\n{len(set(all_not_found))} reactions not found for this set of pathways : \n{', '.join(e for e in set(all_not_found))}")

    # count the number of metabolites whose completion percent is higher than 80% for each strain
    dico_metabo = {}
    metabo_files = os.listdir(output_path + 'metabo_files/')
    for name in metabo_files :
        fr = open(output_path + 'metabo_files/' + name)
        for line in fr :
            line = line[:-1].split('\t')
            if line[0] != 'reaction' :
                filename = line[0]
                if filename not in dico_metabo :
                    dico_metabo[filename] = 0
                percent = line [-1]
                if float(percent) > 80 :
                    dico_metabo[filename] += 1
        fr.close()

    ## completing result files with number of validated pathways and status
    files_created = []
    mkdir(output_path + 'metabo_files_enriched')
    
    for name in metabo_files :  
                   
        with open(output_path + 'metabo_files/' + name, 'r') as fr:
            lines = fr.readlines()
            for i, line in enumerate(lines): 
                listline = line.split("\t")
                if listline[0] != 'reaction' :
                    current_filename = listline[0]
                    if current_filename == "Filename" : 
                        line_extension = f"status\tnb_voies_>_80%"
                    else : 
                        line_extension = f"{filename2status[current_filename]}"
                        try : 
                            line_extension = f"{line_extension}\t{int(dico_metabo[current_filename])}"
                        except :
                            print(f"failed to attribute {current_filename} to nb_voies_>_80%")
                            continue
                    lines[i] = line.replace('\n', f'\t{line_extension}\n')
        new_file = output_path + 'metabo_files_enriched/' + name
        with open(new_file, 'w') as fw:
            for line in lines : 
                fw.write(line)
            files_created.append(name)
    print(f"{len(files_created)} enriched pathway files created :\n{', '.join(e for e in files_created)}")

    #-------------------------------------------------------
        # CREATION OF ASKOMICS FILES
    #-------------------------------------------------------
    if options.asko == True :

        # Writing the tables : Souche, Espece and Genre
        fr1 = open(strain_file)
        mkdir(output_path + 'asko_files')
        fw1 = open(output_path + 'asko_files/' + 'souche.tsv','w')
        fw2 = open(output_path + 'asko_files/' + 'espece.tsv','w')
        fw3 = open(output_path + 'asko_files/' + 'genre.tsv','w')

        fw1.write('Souche\tNom\tassocie@Espece\tStatut\tNombre voies completes a > 80%\n') ## souche.tsv
        fw2.write('Espece\tNom\tassocie@Genre\n')       ## espece.tsv
        fw3.write('Genre\tNom\n')                       ## genre.tsv

        dico_genre = {'Ab':'Anaerobutyricum', 'As':'Anaerostipes','B':'Bifidobacterium','Co':'Coprococcus','C':'Clostridium','Fa':'Faecalibacterium','F':'Fructilactobacillus','E':'Enterococcus','Lb':'Lactobacillus','Lco':'Lactococcus','Leu':'Leuconostoc','Pe':'Pediococus','Pr':'Propionobacterium','St':'Streptococcus','W':'Weissella'}
        esp_list = []
        genre_list = []
        filename2status = {}
        for line in fr1 :
            line = line.split('\t') 
            if line[0] != 'Souche' :
                if options.rename :
                    souche = forbidden(rename(line[0]))
                    filename = forbidden(rename(line[2][:-1]))
                else :
                    souche = forbidden(line[0])
                    filename = forbidden(line[2][:-1])

                if filename in dico_metabo.keys() :
                    statut = line[1]
                    filename2status[filename] = statut      # useful for tsv results' completion
                    espece = souche.split("_")[1]
                    genre = souche.split("_")[0]
                    # Distinction of the two different 'lactis' species
                    if espece == 'lactis' and genre in ['Lco','Lactococcus'] :
                        espece = espece + "_lco"
                    elif espece == 'lactis' and genre in ['Leu','Leuconostoc'] :
                        espece = espece + '_leu'
                    if genre in dico_genre :
                        genre = dico_genre[genre]
                    fw1.write(souche + '\t' + filename + '\t' + espece + '\t' + statut + '\t' + str(dico_metabo[filename]) + '\n')
                    if espece not in esp_list :
                        fw2.write(espece + '\t' + espece + '\t' + genre + '\n')
                        esp_list.append(espece)
                    if genre not in genre_list :
                        fw3.write(genre + '\t' + genre + '\n')
                        genre_list.append(genre)

        fr1.close()
        fw1.close()
        fw2.close()
        fw3.close()
        
        # Modification of the results tables
        # Calculation of the percentages of occurrence of the metabolite pathways, complete at 100% or higher than 80%
        pwy_dict100 = {}
        pwy_dict80 = {}
        list_metabo = []

        ## Writing results for each metabolite file 
        for name in sorted(metabo_files) :
            df = pd.read_csv(output_path + 'metabo_files/' + name ,sep = '\t')
            p,n = df.shape
            df["reaction"] = df["reaction"].map(filename2strain)
            df.rename({'reaction':'associe@Souche'},axis='columns',inplace = True)
            df.insert(0,'Result_'+name[:-4],['result_' + name[:-4] + str(i) for i in range(p)])
            df.insert(1,'associe@Metabolite',[name[:-4] for i in range(p)])

            df.to_csv(output_path + 'asko_files/' + 'result_' + name,sep = '\t', index = False)

            df_percent = df['Completion percent'].values
            df_percent = list(df_percent)
            pwy_dict100[name] = 0
            pwy_dict80[name] = 0
            for i in df_percent[1:] :
                percent = float(i)
                if percent > 80 :
                    pwy_dict80[name] += 1
                    if percent == 100 :
                        pwy_dict100[name] += 1
            list_metabo.append(name)

        ## writing metabolite.tsv    
        fm = open(output_path + 'asko_files/' + 'Metabolite.tsv','w')
        fm.write('Metabolite\tNom\tOccurrence >80%\tOccurrence 100%\n')
        for name in list_metabo :
            fm.write(name[:-4] + '\t' + name[:-4] + '\t' + str(round(pwy_dict80[name]*100/len(df_percent),2)) + '\t' + str(round(pwy_dict100[name]*100/len(df_percent),2)) + '\n')
        fm.close()

if __name__ == "__main__":
    main()
