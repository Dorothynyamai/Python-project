import os
file_list = []
def filename_list():
    user_file = raw_input('please enter a pdb file name: ')
    if user_file in os.listdir('.') and user_file.endswith('.pdb'):

        if user_file not in file_list:
            file_list.append(user_file)

        else:
            print "This file has been read"

    else:
        print "This file is not available"

    new_files = raw_input("Do You want to add a new file ? Yes or No:")
    while new_files.title() == "Yes":
        user_file = raw_input('please enter a file to add to read_list')
        if user_file in os.listdir('.') and user_file.endswith('.pdb'):
            file_list.append(new_file)
        else:
            print "This file is not available"
            new_files = raw_input("Do You want to add a new file ? Yes or No")


    return file_list
#filename_list()
#print read_list


user_option = """select one of the following:
R - read in pdb files
S - search option
W - write option
I - Information option
H - help option
Q - Quit
"""

user_selection = raw_input(user_option)

while user_selection.upper() != "Q":
    if user_selection.upper() == 'R':
        filename_list() #calling the function with list of files. The program checks if the files are in the system before adding them to the list
        print file_list #the list of fies from the file list function
        read_list = []
        file_1 = raw_input("please enter a file from the file_list: ")
        if file_1 in file_list:
            new_file = open(file_1, 'r')
            lines = new_file.readlines()
            new_file.close()
            if lines in read_list:
                print 'The file has been read in'
            else:
                for line in lines:
                    line = line.split()
                    read_list.append(line)
                #print read_list
        else:
            print 'The file you entered is not in the list'


    elif user_selection.upper() == 'W':
        while len(file_list) == 0:
            filename_list()
        else:
            print file_list

        write_out_option = """select one of the following:
        WC - write out coordinate file
        SQ -  write out sequence (FASTA format)
        : """

        user_w_option = raw_input(write_out_option)

        if user_w_option.upper() == 'WC':

            write_sub_option = """select one of the following:
            AA - all atoms
            BA - backbone atoms
            AC - alpha carbon atoms
            : """
            user_sub_option = raw_input(write_sub_option)

            user_file = raw_input('Please enter a file in read_list: ')

            #if user_file in read_list:
            
            #All atom

            if user_sub_option.upper() == 'AA':
                new_file = open(user_file, 'r')
                lines = new_file.readlines()
                new_file.close()
                coor_list = []
                user_file = raw_input('Please enter a file to write the all atom sequence: ')
                for line in lines:
                    if line.startswith('ATOM') or line.startswith("TER"):
                        line = line.strip()
                    coor_list.append(line)

                with open(user_file, 'w') as w:
                    for line in coor_list:
                        w.write(line)
                        w.write('\n')
            
            #Backbone atoms only

            elif user_sub_option.upper() == 'BA':

                ba_file = open(user_file, 'r')
                lines = ba_file.readlines()
                ba_file.close()
                coor_list = []
                user_file = raw_input('Please enter a file to write backbone sequence: ')
                for line in lines:
                    if line.startswith('ATOM') or line.startswith("TER"):
                        line = line.strip()
                        if line[13:15] == 'CA' or line[13:15] == 'C ' or line[13:15] == 'N ' or line[13:15] == 'O ':

                            coor_list.append(line)
                with open(user_file, 'w') as wb:
                    for line in coor_list:
                        wb.write(line)
                        wb.write('\n')
                        
            #Alpha carbon atoms only
            elif user_sub_option.upper() == 'AC':
                c_file = open(user_file, 'r')
                lines = c_file.readlines()
                c_file.close()
                coor_list = []
                user_file = raw_input('Please enter a file to write alpha carbon sequence: ')

                for line in lines:
                    if line.startswith("ATOM") or line.startswith('TER'):
                        line = line.strip()
                        if line[13:15] == 'CA':
                            coor_list.append(line)

                with open(user_file, 'w') as wc:
                    for line in coor_list:
                        wc.write(line)
                        wc.write('\n')

        #Write out sequence (FASTA format)
        
        #SEQRES Sequence
        
        elif user_w_option.upper() == 'SQ':

            write_out_seq_option = """select one of the following:
                s_seq - SEQRES Sequence
                coord_seq - Coordinate Sequence
                al_seq = Alignment Sequence
                : """
            user_seq_option = raw_input(write_out_seq_option)

            user_file = raw_input('Please enter a file in file_list: ')
            seq_list = []
            aa_dict = {'CYS': "C", 'GLN': "Q", 'ILE': "I", 'SER': "S", 'VAL': "V", 'GLY': "G", 'ASN': "N",
                        'PRO': "P", 'LYS': "K", 'THR': "T",'PHE': "F", 'ALA': "A", 'HIS': "H",'MET': "M",
                        'ASP': "D", 'LEU': "L", 'ARG': "R",'TRP': "W", 'GLU': "E", 'TYR': "Y"}


            if user_seq_option.lower() == 's_seq':
                new_file = open(user_file, 'r')
                lines = new_file.readlines()
                new_file.close()
                seqres_file = raw_input('please enter a pdb file name to write the seqres sequence: ')
                #seq_res_list = []
                seq_list = []
                for line in lines:
                    if line.startswith('SEQRES'):
                        line = line.strip()
                        seq_list += (line[19:].split(' '))

                seqRes_seq = ""
                for a in seq_list:
                    seqRes_seq += aa_dict[a]
                with open(seqres_file, 'w') as sq:
                    for line in seqRes_seq:
                        sq.write(line)
                        
                 #Coordinate Sequence       

            elif user_seq_option.lower() == 'coord_seq':
                coord_file = raw_input('please enter a pdb file name to write the coordinate sequence: ')
                coord_seq_list = []
                new_file = open(user_file, 'r')
                lines = new_file.readlines()
                new_file.close()
                atom_pos = ''

                for line in lines:
                    if line.startswith("ATOM") or line.startswith("TER"):
                        line = line.strip()
                        if line[23:26] != atom_pos:
                            coord_seq_list.append(line[17:20])
                            atom_pos = line[23:26]

                coord_list = ''
                for i in coord_seq_list:
                    coord_list += aa_dict[i]

                with open(coord_file, 'w') as c_seq:
                    for line in coord_list:
                        c_seq.write(line)

            #Alignment Sequence
            
            elif user_seq_option.lower() == 'al_seq':
                align_file = raw_input('please enter a pdb file name to write the alignment sequence: ')
                new_file = open(user_file, 'r')
                lines = new_file.readlines()
                new_file.close()
                #to get the sequence list in fasta format
                seq_list = []
                for line in lines:
                    if line.startswith('SEQRES'):
                        line = line.strip()
                        seq_list += (line[19:].split(' '))

                seqRes_seq = ""
                for a in seq_list:
                    seqRes_seq += aa_dict[a]
                #to get the coordinate sequence in fasta format
                coord_seq_list = []
                atom_pos = ''

                for line in lines:
                    if line.startswith("ATOM") or line.startswith("TER"):
                        line = line.strip()
                        if line[23:26] != atom_pos:
                            coord_seq_list.append(line[17:20])
                            atom_pos = line[23:26]

                coord_list = ''
                for i in coord_seq_list:
                    coord_list += aa_dict[i]
                #print coord_list
                #getting aligned sequence
                if len(seqRes_seq) == len(coord_list):
                    print 'No missing amino acids in the sequence'
                elif len(seqRes_seq) != len(coord_list):
                    a = 0
                    while len(seqRes_seq[a]) != len(coord_list[a]):
                        if seqRes_seq[a] != coord_list[a]:
                            coord_list = coord_list[0:a] + "X" + coord_list[a:]
                            a += 1
                        elif seqRes_seq[a] == coord_list[a]:
                            a += 1
                        if a == len(coord_list):
                            break  # means end of the coord_list
                    if len(seqRes_seq) != len(coord_list):
                        aa_diff = len(seqRes_seq) - len(coord_list)
                        coord_list = coord_list + "X" * int(aa_diff)
                        
                        with open(align_file, 'w') as al_seq:
                            for line in coord_list:
                                al_seq.write(line)
                                al_seq.write("\n")

                    elif len(seqRes_seq) == len(coord_list):
                        with open(align_file, 'w') as al_seq:
                            for line in coord_list:
                                al_seq.write(line)
                                al_seq.write("\n")
                                

    elif user_selection.upper() == 'S':
        while len(file_list) == 0:
            filename_list()
        else:
            print file_list
        search_option = """select one of the following:
        mot_search - Motif search
        glyc_search - Glycosylation site search
        : """
        aa_dict = {'CYS': "C", 'GLN': "Q", 'ILE': "I", 'SER': "S", 'VAL': "V", 'GLY': "G", 'ASN': "N",
                    'PRO': "P", 'LYS': "K", 'THR': "T", 'PHE': "F", 'ALA': "A", 'HIS': "H", 'MET': "M",
                    'ASP': "D", 'LEU': "L", 'ARG': "R", 'TRP': "W", 'GLU': "E", 'TYR': "Y"}

        user_sel = raw_input(search_option)

        user_file = raw_input('Please enter a file in file_list: ')
        if user_file in file_list:
            new_file = open(user_file, 'r')
            lines = new_file.readlines()
            new_file.close()
            
            #Motif search
            if user_sel.lower() == 'mot_search':
                user_motif = raw_input('please enter a motif to search: ')

                seq_list = []
                for line in lines:
                    line = line.strip()
                    if line.startswith("SEQRES"):
                        seq_list += (line[19:].split(' '))
                seqres_list = ''
                for i in seq_list:
                    seqres_list += aa_dict[i]

                #index_list = []
                seqres_seq_list = []
                if user_motif.upper() in seqres_list:
                    seqres_seq_list = seqres_list.replace(user_motif.upper(), user_motif.lower())
                    print seqres_seq_list

                    motif_index = []
                    b = 0
                    for a in range(0,len(seqres_seq_list)):
                        if user_motif.lower() == seqres_seq_list[a:a+len(user_motif)]:
                            b += 1
                            motif_index.append(a)
                    print "The number of times the motif appears in the SeQres section is %s." % b
                    print 'The positions of the motif is %s' % motif_index

                else:
                    print "The amino acid motif you entered is not in the sequence"

            #Glycosylation site search
            elif user_sel.lower() == 'glyc_search':
                seq_list = []
                if user_file in file_list:
                    new_file = open(user_file, 'r')
                    lines = new_file.readlines()
                    new_file.close()

                    for line in lines:
                        line = line.strip()
                        if line.startswith("SEQRES"):
                            seq_list += (line[19:].split(' '))
                    seqres_list = ''
                    for i in seq_list:
                        seqres_list += aa_dict[i]

                    import re
                    seqres_seq_list = []
                    glyco_pattern = re.compile('[GYL]A[PFW][TLVGMI]')
                    gly_motif = glyco_pattern.findall(seqres_list)
                    for gly_match in gly_motif:
                        seqres_seq_list = seqres_list.replace(gly_match.upper(), gly_match.lower())


                        gly_index = []
                        b = 0
                        for a in range(0,len(seqres_seq_list)):
                            if gly_match.lower() == seqres_seq_list[a:a+4]:
                                b += 1
                                gly_index.append(a)
                    print seqres_seq_list
                    print "The number of times the glycosylation site motif appears in the SeQres section is %s." % b
                    print 'The index of the glycosilation sites is %s' % gly_index


    elif user_selection.upper() == 'I':
        while len(file_list) == 0:
            filename_list()
        else:
            print file_list
        print_out_info_option = """select one of the following:
        c_seq_info - Display coordinate sequence
        seqres_info - Display SEQRES sequence
        alin_seq_info - Display Alignment sequence
        non_water_lig_info - Display all non-water ligands in the protein
        : """

        aa_dict = {'CYS': "C", 'GLN': "Q", 'ILE': "I", 'SER': "S", 'VAL': "V", 'GLY': "G", 'ASN': "N",
                   'PRO': "P", 'LYS': "K", 'THR': "T", 'PHE': "F", 'ALA': "A", 'HIS': "H", 'MET': "M",
                   'ASP': "D", 'LEU': "L", 'ARG': "R", 'TRP': "W", 'GLU': "E", 'TYR': "Y"}
        user_info_option = raw_input(print_out_info_option)

        if user_info_option == 'c_seq_info':

            user_file = raw_input('Please enter a file in file_list: ')
            for user_file in file_list:
                coord_seq_list = []
                new_file = open(user_file, 'r')
                lines = new_file.readlines()
                new_file.close()
                atom_pos = ''

            for line in lines:
                if line.startswith("ATOM") or line.startswith("TER"):
                    if line[23:26] != atom_pos:
                        coord_seq_list.append(line[17:20])
                        atom_pos = line[23:26]

            coord_list = ''
            for i in coord_seq_list:
                coord_list += aa_dict[i]
            print coord_seq_list

        elif user_info_option == 'seqres_info':
            user_file = raw_input('Please enter a file in file_list: ')
            seq_res_list = []
            new_file = open(user_file, 'r')
            lines = new_file.readlines()
            new_file.close()
            for line in lines:
                seq_list = []
                if line.startswith('SEQRES'):
                    seq_list.append(line.strip())

            for line in seq_list:
                line = line.strip()
                seq_res_list += line[19:].split(' ')

            seqRes_seq = ""
            for a in seq_res_list:
                seqRes_seq += aa_dict[a]

            print seq_list

        elif user_info_option == 'alin_seq_info':
            user_file = raw_input('Please enter a file in file_list: ')

            seq_res_list = []
            new_file = open(user_file, 'r')
            lines = new_file.readlines()
            new_file.close()
            for line in lines:
                seq_list = []
                if line.startswith('SEQRES'):
                    seq_list.append(line.strip())

            for line in seq_list:
                line = line.strip()
                seq_res_list += line[19:].split(' ')

                seqRes_seq = ""
            for a in seq_res_list:
                seqRes_seq += aa_dict[a]

            coord_seq_list = []
            new_file = open(user_file, 'r')
            lines = new_file.readlines()
            new_file.close()
            atom_pos = ''

            for line in lines:
                if line.startswith("ATOM") or line.startswith("TER"):
                    if line[23:26] != atom_pos:
                        coord_seq_list.append(line[17:20])
                        atom_pos = line[23:26]

            coord_list = ''
            for i in coord_seq_list:
                coord_list += aa_dict[i]

            if len(seqRes_seq) == len(coord_list):
                print 'No missing amino acids in the sequence'
            elif len(seqRes_seq) != len(coord_list):
                a = 0
                while len(seqRes_seq) != len(coord_list):
                    if seqRes_seq[a] != coord_list[a]:
                        coord_list = coord_list[0:a] + "X" + coord_list[a:]
                        a += 1
                    elif seqRes_seq[a] == coord_list[a]:
                        a += 1
                    if a == len(coord_list):
                        break #means end of the coord_list
                if len(seqRes_seq) != len(coord_list):
                    aa_diff = len(seqRes_seq) - len(coord_list)
                    coord_list = coord_list + "X" * int(aa_diff)
                    print coord_list
                elif len(seqRes_seq) == len(coord_list):
                    print coord_list


        elif user_info_option == 'non_water_lig_info':

            user_file = raw_input('Please enter a file in file_list: ')
            if user_file in file_list:
                new_file = open(user_file, 'r')
                lines = new_file.readlines()
                new_file.close()
            hetatm_list = []
            for line in lines:
                if line.startswith("HETATM"):
                    line = line.strip()
                    if line[17:20] != "HOH":
                        hetatm_list.append(line)
            print hetatm_list
            
    #Lists all options and gives a brief description of each.        
    elif user_selection.upper() == 'H':


        help_option = """select one of the following:
                        R - read in pdb files
                        S - search option
                        W - write option
                        I - Information option
                        H - help option
                        Q - Quit
                        """

        user_option = raw_input(help_option)

        if user_option.upper() == 'R':
            print "The read option requires the user to enter a PBD file that is in the system which is then then read by the program. If the file the user enters has been read on, the program returns a message that the file has been read on. If the file is not in the system a message is desplayed stating the file is not available. All the file names of the files that have been read on are then stored in a list called read list"
            user_option = raw_input(help_option)
        elif user_option.upper() == 'S':
            print 'This section displays files that have been read in and provides the following suboptions'
            search_option = """ select one of the following:
                            ms - motif search
                            gly_search - Glycosylation site search
                            """
            user_selection = raw_input(search_option)
            if user_selection.lower() == 'ms':
                print "This option allows one to enter in a specific motif then the program prints out the sequence with the matching motifs in lower case and the rest of the sequence in upper case"
                user_selection = raw_input(search_option)
            elif user_selection.lower() == 'gly_search':
                print "This option displays the number of glycosylation sites and their positions in the protein sequence"

                user_option = raw_input(help_option)
        elif user_option.upper() == 'W':
            print "This section has the following sub-options"
            write_out_option = """select one of the following:
                            WC - write out coordinate file
                            SQ -  write out sequence (FASTA format)
                            """
            user_sel = raw_input(write_out_option)
            if user_sel.upper() == 'WC':
                print "The section has the sub-options shown below"
                write_sub_option = """select one of the following:
                            AA - all atoms
                            BA - backbone atoms
                            AC - alpha carbon atoms
                            """
                user_sel = raw_input(write_sub_option)
                if user_sel.upper() == 'AA':
                    print "This section writes all the atoms in the coordinate section to a file"
                    user_sel = raw_input(write_sub_option)
                elif user_sel.upper() == 'BA':
                    print "This section writes only the coordinates of the backbone atoms(CA,C,N,O) to a file"
                    user_sel = raw_input(write_sub_option)
                elif user_sel.upper() == 'AC':
                    print "This section writes only the coordinates for all alpha-carbon atoms(CA) to a file"
                    user_sel = raw_input(write_sub_option)
            elif user_sel.upper() == 'SQ':
                print "This section has the following sub_options"

                write_out_seq_option = """select one of the following:
                            s_seq - SEQRES Sequence
                            coord_seq - Coordinate Sequence
                            al_seq = Alignment Sequence
                            """
                user_sel = raw_input(write_out_seq_option)
                if user_sel.lower() == 's_seq':
                    print "This section writes the SEQRES sequence of the PBD file to a file based on the SEQRES section"
                    user_sel = raw_input(write_out_seq_option)
                elif user_sel.lower() == 'coord_seq':
                    print "This section writes the sequence of the PBD file to a file based on the coordinate section"
                    user_sel = raw_input(write_out_seq_option)
                elif user_sel.lower() == 'al_seq':
                    print "This section writes out a sequence equivalent to the SEQRES section with all the missing residues converted to 'X' characters"
                    user_option = raw_input(help_option)

        elif user_option.upper() == 'I':
            print "This section displays files that have been read in and allows one to choose one file from this list to work on. It contains the following sub_options"

            print_out_info_option = """select one of the following:
                        c_seq_info - Display coordinate sequence
                        seqres_info - Display SEQRES sequence
                        alin_seq_info - Display Alignment sequence
                        non_water_lig_info - Display all non-water ligands in the protein
                        """

            user_sel = raw_input(print_out_info_option)

            if user_sel.lower() == 'c_seq_info':
                print "This section displays the sequence of the PBD file to the screen based on coordinate section"
                user_sel = raw_input(print_out_info_option)

            elif user_sel.lower() == 'seqres_info':
                print "This section displays the SEQRES sequence of the PBD file on the screen based on the SEQRES section"
                user_sel = raw_input(print_out_info_option)
            elif user_sel.lower() == 'alin_seq_info':
                print "This section displays a sequence equivalent to the SEQRES section with all the missing residues converted to 'X' characters to the screen"
                user_sel = raw_input(print_out_info_option)
            elif user_sel.lower() == 'non_water_lig_info':
                print "This section displays the all non- water ligands in the protein based on the HETATM section"
                user_option = raw_input(help_option)

        elif user_option.upper() == 'H':
            print "This option displays all the options in the program and gives a brief description in each"
            user_option = raw_input(help_option)
        elif user_option.upper() == 'Q':
            print "This option allows one to exit from the program"
            break
    elif user_selection.upper() == "Q":
        print "Thank You, you now exit the program"
        quit()
