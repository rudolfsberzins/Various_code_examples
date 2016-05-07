### Made by Rudolfs Berzins, ckj464

import Bio.PDB as PDB
import matplotlib.pyplot as plt
from matplotlib import style
import cmath
import os
from matplotlib.colors import LogNorm
from matplotlib import cm

# Style used for all the plots!
style.use('fivethirtyeight')


def single_protein(name_in_quotes):
    '''By entering a protein PDB ID the function generates 2x2 plot that shows the necessary angle relationships for
     that particular protein. It does so by downloading the protein's structure from PDB, acquiring atom vectors if the
     residue has a C-alpha atom and calculating the dihedral angles. Then it plots the angles with necessary
     specifications on one figure - that way generating 4 plots on 1 image which then is saved in the created folder.'''

    path_to_folder = os.getcwd()  # Necessary for making a folder for saving plots

    # Makes a folder for saving plots
    if not os.path.exists('Ramachandran_Plots'):
        os.makedirs("Ramachandran_Plots")

    # Try is used, because this function catches errors that can occur due to flawed PDB structure files.
    # It will be explained in detail further down.
    try:
        # Download pdb
        pdbl = PDB.PDBList()
        pdbl.retrieve_pdb_file(name_in_quotes, pdir='.')

        # Load pdb file
        parser = PDB.PDBParser()
        structure = parser.get_structure(name_in_quotes, "pdb" + name_in_quotes.lower() + ".ent")

        def calc_phi_psi(structure):
            '''Function makes 3 lists of proteins C-alpha, C and N atom vectors. These lists are then used to calculate
             dihedral angles of proteins.'''

            atom_vector_list_Ca = []
            atom_vector_list_N = []
            atom_vector_list_C = []

            # For-loop for acquiring atom vectors, but only for those residues which have a C-alpha atom.
            for chain in structure.get_chains():
                for res in chain:
                    if res.has_id('CA'):
                        for atom in res:
                            if atom.get_name() == 'N':
                                atom_vector_list_N.append(atom.get_vector())
                            elif atom.get_name() == 'CA':
                                atom_vector_list_Ca.append(atom.get_vector())
                            elif atom.get_name() == 'C':
                                atom_vector_list_C.append(atom.get_vector())
                            else:
                                pass

            len_vec = 0

            ### The if statement compares vector list length between C-alpha vector list and two others, if one of them is
            ### shorter than C-alpha, possibly due to an error in the PDB structure, then the length vector which is
            ### required for calculating dihedral angles is set to be C-alpha which is the same length as other vector lists
            if len(atom_vector_list_Ca) > len(atom_vector_list_C) or len(atom_vector_list_Ca) > len(atom_vector_list_N):
                c_ca = len(atom_vector_list_Ca) - len(atom_vector_list_C)
                n_ca = len(atom_vector_list_Ca) - len(atom_vector_list_N)
                if c_ca == n_ca:
                    len_vec = len(atom_vector_list_Ca) - c_ca
            else:
                len_vec = len(atom_vector_list_Ca)

            dihedral_phi = []
            dihedral_psi = []
            # So we don't include first amino acid which has no phi angle and last amino acid which has no psi angle!
            cut_off = range(1, len_vec - 1)

            # Calculation of phi angles!
            for i in cut_off:
                dihedral_phi.append(PDB.calc_dihedral(atom_vector_list_C[i - 1],
                                                      atom_vector_list_N[i],
                                                      atom_vector_list_Ca[i],
                                                      atom_vector_list_C[i]))

            # Calculation of psi angles!
            for i in cut_off:
                dihedral_psi.append(PDB.calc_dihedral(atom_vector_list_N[i],
                                                      atom_vector_list_Ca[i],
                                                      atom_vector_list_C[i],
                                                      atom_vector_list_N[i + 1]))
            return (dihedral_phi, dihedral_psi)

        dihedral_phi_list, dihedral_psi_list = calc_phi_psi(structure)

        fig = plt.figure(figsize=(10, 10))  # Sets figure size

        fig.suptitle(name_in_quotes, fontsize=14)  # Sets title of the plot as PDB ID

        # Positioning the plots on the figure
        ax1 = plt.subplot2grid((2, 2), (0, 0), rowspan=1, colspan=1)
        ax2 = plt.subplot2grid((2, 2), (0, 1), rowspan=1, colspan=1)
        ax3 = plt.subplot2grid((2, 2), (1, 0), rowspan=1, colspan=1)
        ax4 = plt.subplot2grid((2, 2), (1, 1), rowspan=1, colspan=1)

        dih_phi_len = len(dihedral_phi_list)
        dih_psi_len = len(dihedral_psi_list)

        bin_num = 40

        # Plotting the plots and making them variables for later use in making a unified color bar for all plots
        c1 = ax1.hist2d(dihedral_phi_list[0:dih_phi_len - 1], dihedral_phi_list[1:dih_phi_len], bins=bin_num,
                        norm=LogNorm())
        c2 = ax2.hist2d(dihedral_phi_list[0:dih_phi_len - 1], dihedral_psi_list[1:dih_psi_len], bins=bin_num,
                        norm=LogNorm())
        c3 = ax3.hist2d(dihedral_psi_list[0:dih_psi_len - 1], dihedral_phi_list[1:dih_phi_len], bins=bin_num,
                        norm=LogNorm())
        c4 = ax4.hist2d(dihedral_psi_list[0:dih_psi_len - 1], dihedral_psi_list[1:dih_psi_len], bins=bin_num,
                        norm=LogNorm())

        lim = [-cmath.pi, cmath.pi, -cmath.pi, cmath.pi]  # Limits of the axes

        ax1.axis(lim)
        ax2.axis(lim)
        ax3.axis(lim)
        ax4.axis(lim)

        # Unified font for all plots
        font_dict = {'family': 'times new roman',
                     'size': 15}

        # Setting labels for axes
        ax1.set_xlabel('$\phi$(i)', fontdict=font_dict)
        ax1.set_ylabel('$\phi$(i+1)', fontdict=font_dict)
        ax2.set_xlabel('$\phi$(i)', fontdict=font_dict)
        ax2.set_ylabel('$\psi$(i+1)', fontdict=font_dict)
        ax3.set_xlabel('$\psi$(i)', fontdict=font_dict)
        ax3.set_ylabel('$\phi$(i+1)', fontdict=font_dict)
        ax4.set_xlabel('$\psi$(i)', fontdict=font_dict)
        ax4.set_ylabel('$\psi$(i+1)', fontdict=font_dict)

        # Acquiring an array of observation count in a specific bin
        count1 = c1[0]
        count2 = c2[0]
        count3 = c3[0]
        count4 = c4[0]

        # Averaging the count in each plot based on mean maximum value between the plots. This is used for color bar
        count = (count1.max() + count2.max() + count3.max() + count4.max()) / 4

        # Creates a unified color bar for all plots, with min value of 1 observations at a position and a max value of
        # the calculated average max value between the 4 plots.
        fig.subplots_adjust(right=0.8)
        a = [1, count]
        m = cm.ScalarMappable(cmap=cm.jet)
        m.set_array(a)
        cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
        fig.colorbar(m, cax=cbar_ax,
                     label='Average height of bins, i.e. observation count (bin number ' + str(bin_num) + ')')

        # Saves figure with PDB ID as the picture name
        fig.savefig("Ramachandran_Plots" + "\\" + name_in_quotes + ".png")

    ### This section is required because some PDB have error in them.
    ### The next lines of code are catching these errors so they don't hinder the work with a large data set if they
    ### appear. The code makes a .txt file in the same folder as the plot images and write a specific text when each
    ### error occurs and also pointing to the PDB file with which the error occured.
    except (IndexError):
        os.chdir(path_to_folder + '\\Ramachandran_Plots')
        flawed_files = open('Flawed_PDB_files.txt', 'a')
        flawed_files.write(
                '' + name_in_quotes + ' because a Index error occurred while using this PDB file. Please check the file! \n')
        os.chdir(path_to_folder)
        pass
    except (TypeError):
        os.chdir(path_to_folder + '\\Ramachandran_Plots')
        flawed_files = open('Flawed_PDB_files.txt', 'a')
        flawed_files.write(
                '' + name_in_quotes + ' because a Type error occurred while using this PDB file. Please check the file! \n')
        os.chdir(path_to_folder)
        pass
    except (NameError):
        os.chdir(path_to_folder + '\\Ramachandran_Plots')
        flawed_files = open('Flawed_PDB_files.txt', 'a')
        flawed_files.write(
                '' + name_in_quotes + ' because a Name error occurred while using this PDB file. Please check the file! \n')
        os.chdir(path_to_folder)
        pass
    except (ValueError):
        os.chdir(path_to_folder + '\\Ramachandran_Plots')
        flawed_files = open('Flawed_PDB_files.txt', 'a')
        flawed_files.write(
                '' + name_in_quotes + ' because a Value error occurred while using this PDB file. Please check the file! \n')
        os.chdir(path_to_folder)
        pass
    except (PermissionError):
        os.chdir(path_to_folder + '\\Ramachandran_Plots')
        flawed_files = open('Flawed_PDB_files.txt', 'a')
        flawed_files.write(
                '' + name_in_quotes + ' because a Permission error occurred while using this PDB file. Please check the file! \n')
        os.chdir(path_to_folder)
        pass


def list_of_structures(path_to_folder, foldername):
    '''By entering a path to the folder containing PDB structures and the actual folder name, the function makes a list of
     all the protein angles combined and makes a Full Ramachandran plot from all the protein angles in the folder,
     based on the specific plotting requirements. It does this by acquiring atom vectors if the residue has a C-alpha atom
     and calculating the dihedral angles. Then it plots the proteins angles with necessary specifications on one figure -
     that way generating 4 plots on 1 image which then is saved in the created folder. After that the function generates
     a traditional Ramachandran plot (Phi/Psi). The function can also generates 2x2 plot that shows the necessary
     angle relationships for each protein structure in the folder. You have to uncomment a section of a code to enable it.'''

    # Changes the working directory for the one where the folder with structures is.
    os.chdir(path_to_folder)
    list_of_proteins = []

    # Makes a folder for saving plots
    if not os.path.exists('Ramachandran_Plots'):
        os.makedirs("Ramachandran_Plots")

    # Takes out individual protein names from the folder of structures
    for name in os.listdir(foldername):
        list_of_proteins.append(name)

    def full_ramachandran(list_of_prot):
        '''Function makes two list that will will be appended with dihedral values of individual proteins. Protein
        structures will be taken from the structure folder. As an argument it takes the list of proteins made previously'''
        full_dih_phi_list = []
        full_dih_psi_list = []
        path_to_folder = os.getcwd()

        # This works similar to the single protein function, but instead of a separate function we use a for-loop.
        for protein in list_of_prot:
            # Try is used, because this function catches errors that can occur due to flawed PDB structure files.
            # It will be explained in detail further down.
            try:
                os.chdir(path_to_folder + '\\' + foldername)

                # Load pdb file
                parser = PDB.PDBParser()
                structure = parser.get_structure(protein, protein.lower())

                def calc_phi_psi(structure):
                    '''Function makes 3 lists of proteins C-alpha, C and N atom vectors. These lists are then used to calculate
                     dihedral angles of proteins.'''

                    atom_vector_list_Ca = []
                    atom_vector_list_N = []
                    atom_vector_list_C = []

                    # For-loop for acquiring atom vectors, but only for those residues which have a C-alpha atom.
                    for chain in structure.get_chains():
                        for res in chain:
                            if res.has_id('CA'):
                                for atom in res:
                                    if atom.get_name() == 'N':
                                        atom_vector_list_N.append(atom.get_vector())
                                    elif atom.get_name() == 'CA':
                                        atom_vector_list_Ca.append(atom.get_vector())
                                    elif atom.get_name() == 'C':
                                        atom_vector_list_C.append(atom.get_vector())
                                    else:
                                        pass

                    len_vec = 0

                    ### The if statement compares vector list length between C-alpha vector list and two others, if one of them is
                    ### shorter than C-alpha, possibly due to an error in the PDB structure, then the length vector which is
                    ### required for calculating dihedral angles is set to be C-alpha which is the same length as other vector lists
                    if len(atom_vector_list_Ca) > len(atom_vector_list_C) or len(atom_vector_list_Ca) > len(
                            atom_vector_list_N):
                        c_ca = len(atom_vector_list_Ca) - len(atom_vector_list_C)
                        n_ca = len(atom_vector_list_Ca) - len(atom_vector_list_N)
                        if c_ca == n_ca:
                            len_vec = len(atom_vector_list_Ca) - c_ca
                    else:
                        len_vec = len(atom_vector_list_Ca)

                    dihedral_phi = []
                    dihedral_psi = []

                    # So we don't include first amino acid which has no phi angle and last amino acid which has no psi angle!
                    cut_off = range(1, len_vec - 1)

                    # Calculation of phi angles!
                    for i in cut_off:
                        dihedral_phi.append(PDB.calc_dihedral(atom_vector_list_C[i - 1],
                                                              atom_vector_list_N[i],
                                                              atom_vector_list_Ca[i],
                                                              atom_vector_list_C[i]))

                    # Calculation of psi angles!
                    for i in cut_off:
                        dihedral_psi.append(PDB.calc_dihedral(atom_vector_list_N[i],
                                                              atom_vector_list_Ca[i],
                                                              atom_vector_list_C[i],
                                                              atom_vector_list_N[i + 1]))
                    return (dihedral_phi, dihedral_psi)

                dihedral_phi_list, dihedral_psi_list = calc_phi_psi(structure)

                # Adds the acquired lists to the two lists mentioned previously
                full_dih_phi_list.extend(dihedral_phi_list)
                full_dih_psi_list.extend(dihedral_psi_list)

                #### UNCOMMENT ALL OF THE COMMENTED SECTION BELOW TO GET INDIVIDUAL RAMACHANDRAN PLOTS FOR EACH PROTEIN
                #### IN THE DATA SET (FOLDER)!!! NOT RECOMMENDED FOR LARGE SETS OF SEQUENCES, LIKE 1000+

                # fig = plt.figure(figsize=(10, 10))  # Sets figure size
                #
                # fig.suptitle(protein, fontsize=14)  # Sets title of the plot as PDB ID
                #
                # # Positioning the plots on the figure
                # ax1 = plt.subplot2grid((2, 2), (0, 0), rowspan=1, colspan=1)
                # ax2 = plt.subplot2grid((2, 2), (0, 1), rowspan=1, colspan=1)
                # ax3 = plt.subplot2grid((2, 2), (1, 0), rowspan=1, colspan=1)
                # ax4 = plt.subplot2grid((2, 2), (1, 1), rowspan=1, colspan=1)
                #
                # dih_phi_len = len(dihedral_phi_list)
                # dih_psi_len = len(dihedral_psi_list)
                #
                # bin_num = 40
                #
                # # Plotting the plots and making them variables for later use in making a unified color bar for all plots
                # c1 = ax1.hist2d(dihedral_phi_list[0:dih_phi_len - 1], dihedral_phi_list[1:dih_phi_len], bins=bin_num,
                #                 norm=LogNorm())
                # c2 = ax2.hist2d(dihedral_phi_list[0:dih_phi_len - 1], dihedral_psi_list[1:dih_psi_len], bins=bin_num,
                #                 norm=LogNorm())
                # c3 = ax3.hist2d(dihedral_psi_list[0:dih_psi_len - 1], dihedral_phi_list[1:dih_phi_len], bins=bin_num,
                #                 norm=LogNorm())
                # c4 = ax4.hist2d(dihedral_psi_list[0:dih_psi_len - 1], dihedral_psi_list[1:dih_psi_len], bins=bin_num,
                #                 norm=LogNorm())
                #
                # lim = [-cmath.pi, cmath.pi, -cmath.pi, cmath.pi]  # Limits of the axes
                #
                # ax1.axis(lim)
                # ax2.axis(lim)
                # ax3.axis(lim)
                # ax4.axis(lim)
                #
                # # Unified font for all plots
                # font_dict = {'family': 'times new roman',
                #              'size': 15}
                #
                # # Setting labels for axes
                # ax1.set_xlabel('$\phi$(i)', fontdict=font_dict)
                # ax1.set_ylabel('$\phi$(i+1)', fontdict=font_dict)
                # ax2.set_xlabel('$\phi$(i)', fontdict=font_dict)
                # ax2.set_ylabel('$\psi$(i+1)', fontdict=font_dict)
                # ax3.set_xlabel('$\psi$(i)', fontdict=font_dict)
                # ax3.set_ylabel('$\phi$(i+1)', fontdict=font_dict)
                # ax4.set_xlabel('$\psi$(i)', fontdict=font_dict)
                # ax4.set_ylabel('$\psi$(i+1)', fontdict=font_dict)
                #
                # # Acquiring an array of observation count in a specific bin
                # count1 = c1[0]
                # count2 = c2[0]
                # count3 = c3[0]
                # count4 = c4[0]
                #
                # # Averaging the count in each plot based on mean maximum value between the plots. This is used for color bar
                # count = (count1.max() + count2.max() + count3.max() + count4.max()) / 4
                #
                # # Creates a unified color bar for all plots, with min value of 1 observations at a position and a max value of
                # # the calculated average max value between the 4 plots.
                # fig.subplots_adjust(right=0.8)
                # a = [1, count]
                # m = cm.ScalarMappable(cmap=cm.jet)
                # m.set_array(a)
                # cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
                # fig.colorbar(m, cax=cbar_ax,
                #              label='Average height of bins, i.e. observation count (bin number ' + str(bin_num) + ')')
                #
                # # Saves figure with PDB ID as the picture name
                # os.chdir(path_to_folder)
                # fig.savefig("Ramachandran_Plots" + "\\" + protein + ".png")
                # fig.clf()

            ### This section is required because some PDB have error in them.
            ### The next lines of code are catching these errors so they don't hinder the work with a large data set if they
            ### appear. The code makes a .txt file in the same folder as the plot images and write a specific text when each
            ### error occurs and also pointing to the PDB file with which the error occured.
            except (IndexError):
                os.chdir(path_to_folder + '\\Ramachandran_Plots')
                flawed_files = open('Flawed_PDB_files.txt', 'a')
                flawed_files.write(
                        '' + protein + ' because a Index error occurred while using this PDB file. Please check the file! \n')
                os.chdir(path_to_folder)
                pass
            except (TypeError):
                os.chdir(path_to_folder + '\\Ramachandran_Plots')
                flawed_files = open('Flawed_PDB_files.txt', 'a')
                flawed_files.write(
                        '' + protein + ' because a Type error occurred while using this PDB file. Please check the file! \n')
                os.chdir(path_to_folder)
                pass
            except (NameError):
                os.chdir(path_to_folder + '\\Ramachandran_Plots')
                flawed_files = open('Flawed_PDB_files.txt', 'a')
                flawed_files.write(
                        '' + protein + ' because a Name error occurred while using this PDB file. Please check the file! \n')
                os.chdir(path_to_folder)
                pass
            except (ValueError):
                os.chdir(path_to_folder + '\\Ramachandran_Plots')
                flawed_files = open('Flawed_PDB_files.txt', 'a')
                flawed_files.write(
                        '' + protein + ' because a Value error occurred while using this PDB file. Please check the file! \n')
                os.chdir(path_to_folder)
                pass
            except (PermissionError):
                os.chdir(path_to_folder + '\\Ramachandran_Plots')
                flawed_files = open('Flawed_PDB_files.txt', 'a')
                flawed_files.write(
                        '' + protein + ' because a Permission error occurred while using this PDB file. Please check the file! \n')
                os.chdir(path_to_folder)
                pass

        return (full_dih_phi_list, full_dih_psi_list)

    full_dih_phi_list, full_dih_psi_list = full_ramachandran(list_of_proteins)

    fig = plt.figure(figsize=(10, 10))  # Sets figure size

    fig.suptitle("Full Ramachandran Plots", fontsize=16)  # Sets title

    # Positioning the plots on the figure
    ax1 = plt.subplot2grid((2, 2), (0, 0), rowspan=1, colspan=1)
    ax2 = plt.subplot2grid((2, 2), (0, 1), rowspan=1, colspan=1)
    ax3 = plt.subplot2grid((2, 2), (1, 0), rowspan=1, colspan=1)
    ax4 = plt.subplot2grid((2, 2), (1, 1), rowspan=1, colspan=1)

    dih_phi_len = len(full_dih_phi_list)
    dih_psi_len = len(full_dih_psi_list)

    bin_num = 360

    # Plotting the plots and making them variables for later use in making a unified color bar for all plots
    c1 = ax1.hist2d(full_dih_phi_list[0:dih_phi_len - 1], full_dih_phi_list[1:dih_phi_len], bins=bin_num,
                    norm=LogNorm())
    c2 = ax2.hist2d(full_dih_phi_list[0:dih_phi_len - 1], full_dih_psi_list[1:dih_psi_len], bins=bin_num,
                    norm=LogNorm())
    c3 = ax3.hist2d(full_dih_psi_list[0:dih_psi_len - 1], full_dih_phi_list[1:dih_phi_len], bins=bin_num,
                    norm=LogNorm())
    c4 = ax4.hist2d(full_dih_psi_list[0:dih_psi_len - 1], full_dih_psi_list[1:dih_psi_len], bins=bin_num,
                    norm=LogNorm())

    lim = [-cmath.pi, cmath.pi, -cmath.pi, cmath.pi]  # Limits of the axes

    ax1.axis(lim)
    ax2.axis(lim)
    ax3.axis(lim)
    ax4.axis(lim)

    # Unified font for all plots
    font_dict = {'family': 'times new roman',
                 'size': 15}

    # Setting labels for axes
    ax1.set_xlabel('$\phi$(i)', fontdict=font_dict)
    ax1.set_ylabel('$\phi$(i+1)', fontdict=font_dict)
    ax2.set_xlabel('$\phi$(i)', fontdict=font_dict)
    ax2.set_ylabel('$\psi$(i+1)', fontdict=font_dict)
    ax3.set_xlabel('$\psi$(i)', fontdict=font_dict)
    ax3.set_ylabel('$\phi$(i+1)', fontdict=font_dict)
    ax4.set_xlabel('$\psi$(i)', fontdict=font_dict)
    ax4.set_ylabel('$\psi$(i+1)', fontdict=font_dict)

    # Acquiring an array of observation count in a specific bin
    count1 = c1[0]
    count2 = c2[0]
    count3 = c3[0]
    count4 = c4[0]

    # Averaging the count in each plot based on mean maximum value between the plots. This is used for color bar
    count = (count1.max() + count2.max() + count3.max() + count4.max()) / 4

    # Creates a unified color bar for all plots, with min value of 1 observations at a position and a max value of
    # the calculated average max value between the 4 plots.
    fig.subplots_adjust(right=0.8)
    a = [1, count]
    m = cm.ScalarMappable(cmap=cm.jet)
    m.set_array(a)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.colorbar(m, cax=cbar_ax,
                 label='Average height of bins, i.e. observation count (bin number ' + str(bin_num) + ')')

    # Saves figure
    os.chdir(path_to_folder)
    fig.savefig('Ramachandran_Plots\Ramachandran_Full_' + foldername + '.png')
    fig.clf()

    # Creates a regular Ramachandran plot (phi/psi), plots it and saves it.
    c5 = plt.hist2d(full_dih_phi_list, full_dih_psi_list, bins=bin_num, norm=LogNorm())
    count5 = c5[0]
    a = [1, count5.max()]
    m = cm.ScalarMappable(cmap=cm.jet)
    m.set_array(a)
    plt.title('Correct Ramachandran plot - $\phi/\psi$')
    plt.xlabel('$\phi$')
    plt.ylabel('$\psi$')
    plt.colorbar(m, label='Height of bins, i.e. observation count (bin number ' + str(bin_num) + ')')
    plt.savefig('Ramachandran_Plots\Ramachandran_Correct_' + foldername + '.png')


# list_of_structures('path to the folder with structures', 'name of the folder')
single_protein("1ENH")
