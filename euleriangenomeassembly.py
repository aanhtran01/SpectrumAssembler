# -*- coding: utf-8 -*-

import argparse
from copy import deepcopy

'''
given a spectrum for a k that is smaller than the reads, a de Bruijn graph is generated to reconstruct a reference genome 

the output is a file contains the headers of the reads sorted in the order that they appear in the genome
'''

#set up a command line interface with argparse, create argparse objects
parser = argparse.ArgumentParser(description='Process spectrum.')
parser.add_argument('spectrum_file', type=str, help='path to spectrum file')


args = parser.parse_args()

spectrum_file = args.spectrum_file

#function to read in fasta files 
def read_fasta_file(filename):
    with open(filename, 'r') as f:
        read_id = None
        read_seq = ""
        for line in f:
            if line.startswith(">"):
                if read_id is not None:
                    yield (read_id, read_seq)
                read_id = line.strip()[1:]
                read_seq = ""
            else:
                read_seq += line.strip()
        if read_id is not None:
            yield (read_id, read_seq)

class Euler:
    def __init__(self, adjacency_list):
        """
        Initializes the Euler class with an adjacency list.

        Parameters:
        - adjacency_list: A dictionary representing the adjacency list of the de Bruijn graph.
        """
        self.adjacency_list = adjacency_list  # Assign the adjacency list
        self.num_nodes = len(self.adjacency_list)  # Calculate the number of nodes in the graph
        self.current_position = {}  # Tracks the current position in the exploration of each node
        self.path = []  # Stores the Eulerian path
        self.unbalanced_nodes = []  # Stores the unbalanced nodes in the graph
        self.nodes_with_unused_edges = {}  # Stores nodes that have unused edges in the current path
        self.num_unexplored_edges = 0  # Counts the number of unexplored edges in the graph
        self.in_degree = {}  # Stores the in-degree of each node
        self.out_degree = {}  # Stores the out-degree of each node
        self.update_adjacency_list()  # Update the adjacency list and calculate degrees

    def update_adjacency_list(self):
    # Update the in-degree and out-degree for each node in the adjacency list
        for node, neighbors in self.adjacency_list.items():
            self.in_degree[node] = self.in_degree.get(node, 0)  # Get the current in-degree of the node (default to 0 if not present)
            for neighbor in neighbors:
                self.in_degree[neighbor] = 1 + self.in_degree.get(neighbor, 0)  # Increment the in-degree of the neighbor node
            self.current_position[node] = 0  # Set the current position of the node to 0
            num_neighbors = len(neighbors)  # Get the number of neighbors of the node
            self.out_degree[node] = num_neighbors  # Set the out-degree of the node to the number of neighbors
            self.num_unexplored_edges = self.num_unexplored_edges + num_neighbors  # Increment the count of unexplored edges

    def add_missing_edge(self):
    # Add a missing edge between unbalanced nodes to make the graph balanced
        if isinstance(self.adjacency_list, dict):  # Check if the adjacency_list is a dictionary
            for node in self.adjacency_list.keys():
                if self.in_degree[node] != self.out_degree[node]:
                    if self.in_degree[node] < self.out_degree[node]:
                        self.unbalanced_nodes.append(node)  # Add the unbalanced node to the end of the list
                    else:
                        self.unbalanced_nodes.insert(0, node)  # Add the unbalanced node to the beginning of the list
            if self.unbalanced_nodes:
                self.update_adjacency_list_for_unbalanced_nodes()  # Update the adjacency list for unbalanced nodes
            return
        for node in range(self.num_nodes):  # Iterate over the range of node indices
            if self.in_degree[node] != self.out_degree[node]:
                if self.in_degree[node] < self.out_degree[node]:
                    self.unbalanced_nodes.append(node)  # Add the unbalanced node to the end of the list
                else:
                    self.unbalanced_nodes.insert(0, node)  # Add the unbalanced node to the beginning of the list
        if self.unbalanced_nodes:
            self.update_adjacency_list_for_unbalanced_nodes()  # Update the adjacency list for unbalanced nodes


    def update_adjacency_list_for_unbalanced_nodes(self):
        """
        Updates the adjacency list for unbalanced nodes by adding a missing edge.

        This method increases the out-degree of the first unbalanced node and the in-degree of the second unbalanced node by 1.
        It also adds the second unbalanced node to the adjacency list of the first unbalanced node.
        """
        self.out_degree[self.unbalanced_nodes[0]] = self.out_degree[self.unbalanced_nodes[0]] + 1  # Increase the out-degree of the first unbalanced node by 1
        self.in_degree[self.unbalanced_nodes[1]] = self.in_degree[self.unbalanced_nodes[1]] + 1  # Increase the in-degree of the second unbalanced node by 1
        self.adjacency_list[self.unbalanced_nodes[0]].append(self.unbalanced_nodes[1])  # Add the second unbalanced node to the adjacency list of the first unbalanced node

        

    def explore_node(self, start_node):
        """
        Explores a node and its neighbors in the graph.

        This method traverses the graph from the given start_node, adding nodes to the path and updating the current position.
        It keeps track of nodes with unused edges and decreases the count of unexplored edges in the graph.
        """
        self.path.append(start_node)  # Add the start_node to the path
        current_max_pos = self.out_degree[start_node]  # Get the maximum position for the start_node
        current_pos = self.current_position[start_node]  # Get the current position of the start_node

        while current_pos < current_max_pos:
            self.current_position[start_node] = current_pos + 1  # Update the current position of the start_node

            if current_pos + 1 < current_max_pos:
                self.nodes_with_unused_edges[start_node] = len(self.path) - 1  # Record the position of the start_node with unused edges
            else:
                if start_node in self.nodes_with_unused_edges:
                    del self.nodes_with_unused_edges[start_node]  # Remove the start_node from nodes with unused edges

            neighbor = self.adjacency_list[start_node][current_pos]  # Get the next neighbor of the start_node
            self.path.append(neighbor)  # Add the neighbor to the path
            start_node = neighbor  # Update the start_node for the next iteration
            current_pos = self.current_position[start_node]  # Update the current position of the new start_node
            current_max_pos = self.out_degree[start_node]  # Update the maximum position for the new start_node
            self.num_unexplored_edges = self.num_unexplored_edges - 1  # Decrease the count of unexplored edges in the graph

    def update_path(self, start_pos):
        """
        Updates the path based on the start position.

        This method rearranges the path by shifting the elements starting from the given start position.
        It also updates the positions of nodes with unused edges in the updated path.
        """
        path_length = len(self.path) - 1  # Calculate the length of the path
        self.path = self.path[start_pos:path_length] + self.path[:start_pos]  # Rearrange the path based on the start position

        for node, pos in self.nodes_with_unused_edges.items():
            # Update the positions of nodes with unused edges in the updated path
            self.nodes_with_unused_edges[node] = (pos + path_length - start_pos) if pos < start_pos else (pos - start_pos)

    
    def eulerian_path(self):
        """
        Finds the Eulerian path in the de Bruijn graph.

        This method constructs the Eulerian path by exploring nodes and updating the path and positions.
        It handles both balanced and unbalanced graphs.
        """
        self.add_missing_edge()  # Add missing edges to balance the graph

        if isinstance(self.adjacency_list, dict):
            node, neighbor_list = self.adjacency_list.popitem()
            self.adjacency_list[node] = neighbor_list
            self.explore_node(node)  # Explore the first node in the adjacency list
        else:
            self.explore_node(0)  # Explore the first node in the range of node indices

        while self.num_unexplored_edges > 0:
            node, pos = self.nodes_with_unused_edges.popitem()
            self.update_path(pos)  # Update the path based on the position
            self.explore_node(node)  # Explore the node with unused edges

        if self.unbalanced_nodes:
            # Handle unbalanced nodes by finding the position where the unbalanced nodes appear consecutively in the path
            for i in range(len(self.path) - 1):
                if self.path[i] == self.unbalanced_nodes[0] and self.path[i + 1] == self.unbalanced_nodes[1]:
                    self.update_path(i + 1)  # Update the path based on the position
                    break

        return self.path  # Return the Eulerian path



class GenomeAssembly:
  
    def __init__(self, filename):
        """
        Initializes the GenomeAssembly object by reading the de Bruijn graph from a file, finding the Eulerian path,
        and reconstructing the genome from the path.
        """
        self.adjacency_list = self.read_de_bruijn(filename)  # Read the de Bruijn graph from the file
        self.path = Euler(self.adjacency_list).eulerian_path()  # Find the Eulerian path in the de Bruijn graph
        self.reconstructed_genome = self.reconstruct_from_path(self.path)  # Reconstruct the genome from the path

    def reconstruct_from_path(self, path):
        """
        Reconstructs the genome sequence from the given path.

        This method concatenates the first character of each sequence in the path, except for the first sequence.
        """
        return path[0] + ''.join(seq[-1] for seq in path[1:])

    def save_result(self):
        """
        Saves the reconstructed genome as the result.

        This method returns the reconstructed genome sequence.
        """
        return self.reconstructed_genome

    def __str__(self):
        """
        Returns the string representation of the reconstructed genome.

        This method allows the GenomeAssembly object to be printed as a string.
        """
        return self.reconstructed_genome

    def __iadd__(self, other):
        """
        Implements the += operator for GenomeAssembly objects.

        This method appends the given string to the reconstructed genome.
        """
        if isinstance(other, str):
            self.reconstructed_genome += other
            return self
        else:
            raise TypeError("Unsupported operand type for +=: 'GenomeAssembly' and '{}'".format(type(other).__name__))

    def __len__(self):
        """
        Returns the length of the reconstructed genome sequence.

        This method allows the len() function to be used on a GenomeAssembly object, returning the length of the genome.
        """
        return len(self.reconstructed_genome)

    def read_de_bruijn(self, filename):
        """
        Reads a FASTA file containing sequences and constructs the de Bruijn graph.

        This method reads the sequences from the FASTA file and constructs the de Bruijn graph by creating adjacency lists
        or each k-1 length prefix and appending the corresponding suffixes. The adjacency lists are stored in a dictionary,
        where the keys are prefixes and the values are lists of suffixes.

        Returns:
        - adjacency_db: A dictionary representing the de Bruijn graph adjacency list.
        """
        adjacency_db = {}  # Dictionary to store the adjacency lists of the de Bruijn graph
        for read_id, read_seq in read_fasta_file(spectrum_file):
            k = len(read_seq) - 1
            prefix = read_seq[:k]  # Extract the k-1 length prefix
            suffix = read_seq[1:]  # Extract the corresponding suffix
            adjacency_db.setdefault(prefix, []).append(suffix)  # Append the suffix to the adjacency list of the prefix
            if suffix not in adjacency_db:
                adjacency_db[suffix] = []  # Create an empty adjacency list for the suffix if it doesn't exist
        return adjacency_db



#BurrowsWheeler  Construction

class suffix_array:
    def __init__(self, text):
        # Constructor that builds the suffix array for a given text string
        self.suff_arr = self._make_suffix_array(text)


    def _sort_array(self, S):
        # Method that sorts the characters of the string S
        l = len(S)
        arrangement = [0] * l
        num = {}
        
        # Count the occurrences of each character in S
        for i in range(l):
            num[S[i]] = num.get(S[i], 0) + 1
        
        # Sort the characters in ascending arrangement
        char_list = sorted(num.keys())
        
        # Compute the starting position of each character in the sorted arrangement
        prev_char = char_list[0]
        for char in char_list[1:]:
            num[char] =  num[char] + num[prev_char]
            prev_char = char
        
        # Compute the arrangement of each suffix based on the starting character
        for i in range(l-1, -1, -1):
            c = S[i]
            num[c] = num[c] - 1
            arrangement[num[c]] = i
        
        return arrangement

    def class_character_arrangement(self, S, arrangement):
        # Method that computes the character classes for each position in the suffix array
        l = len(S)
        class_characters = [0] * l
        class_characters[arrangement[0]] = 0
        
        # Assign the class to each suffix based on whether the starting character is the same as the previous suffix
        for i in range(1, l):
            if S[arrangement[i]] != S[arrangement[i-1]]:
                class_characters[arrangement[i]] = class_characters[arrangement[i-1]] + 1
            else:
                class_characters[arrangement[i]] = class_characters[arrangement[i-1]]
        
        return class_characters

    def _double_suffix_sort(self, S, L, arrangement, class_characters):
        # Method that sorts the doubled suffixes
        string_length = len(S)
        num = [0] * string_length
        new_arrangement = [0] * string_length
        
        # Count the occurrences of each class in the first half of the doubled suffixes
        for i in range(string_length):
            num[class_characters[i]] = num[class_characters[i]] + 1
        
        # Compute the starting position of each class in the sorted arrangement
        for j in range(1, string_length):
            num[j] = num[j] + num[j-1]
        
        # Sort the doubled suffixes based on their second half
        for i in range(string_length-1, -1, -1):
            start = (arrangement[i]-L+string_length) % string_length
            cl = class_characters[start]
            num[cl] = num[cl] - 1
            new_arrangement[num[cl]] = start
        
        return new_arrangement
    
    def _update_classes(self, new_arrangement, class_characters, L):
      # This function updates the character classes based on a new arrangement of indices.
      # new_arrangement: the new arrangement of indices for the string
      # class_characters: a list containing the character classes of the string
      # L: the length of substrings used for sorting
    
      n = len(new_arrangement)
     # n is the length of the new arrangement
    
      class_new = [0] * n
      # Create a new list of character classes with n elements initialized to 0
    
      class_new[new_arrangement[0]] = 0
      # The character class of the first element in the new arrangement is always 0
    
      for i in range(1, n):
          prev = new_arrangement[i-1]
          curr = new_arrangement[i]
          mid = curr + L
          mid_prev = (prev + L) % n
          # Define curr, prev, mid, and mid_prev for easier readability
        
          # Compare the character classes of two adjacent elements in the new arrangement
          # and two adjacent substrings starting at those elements, respectively.
          # If they're different, increment the character class.
          if class_characters[curr] != class_characters[prev] or class_characters[mid] != class_characters[mid_prev]:
              class_new[curr] = class_new[prev] + 1
          else:
              class_new[curr] = class_new[prev]
    
      # Return the updated character classes
      return class_new
    
    def _make_suffix_array(self, S):
      # This function builds the suffix array for a given string S
      # S: the input string
    
      string_length = len(S)
      # The length of S
    
      arrangement = self._sort_array(S)
      # Sort the characters of S and store the resulting arrangement
    
      class_characters = self.class_character_arrangement(S, arrangement)
      # Compute the character classes of S and store them
    
      L = 1
      while L < string_length:
          arrangement = self._double_suffix_sort(S, L, arrangement, class_characters)
          class_characters = self._update_classes(arrangement, class_characters, L)
          L = 2 * L
      # Repeat the process with longer substrings until L >= string_length.
    
      # Return the final suffix array
      return arrangement


    def get_suffix_array(self):
      # This function returns the suffix array of the input string.
    
      return self.suff_arr
      # Return the suffix array stored in self.sa


class BurrowsWheeler:
    def __init__(self, reference_genome):
        self.burrows_wheeler = self.burrows_wheelerFromsuffix_array(reference_genome)
    # Initialize the BurrowsWheeler by calling burrows_wheelerFromsuffix_array with a reference genome
    
    def burrows_wheelerTransform(self, text):
        # This function calculates the Burrows-Wheeler Transform
        # text: the input string
        
        # Generate all possible transfrom of the input string
        transfrom = [text[i:]+text[:i] for i in range(len(text))]
        # Sort the transfrom lexicographically and concatenate their last characters
        burrows_wheeler = ''.join([m[-1] for m in sorted(transfrom)])
        
        return burrows_wheeler

    def burrows_wheelerFromsuffix_array(self, text):
       # This function calculates the Burrows-Wheeler Transform using suffix arrays.
        suff_arr = suffix_array(text).get_suffix_array()
        return ''.join([text[(suff_arr[i]+len(text)-1)%len(text)] for i in range(len(text))])


def HammingDistance(seq1, seq2):
    return len([i for i in range(len(seq1)) if seq1[i] != seq2[i]])


#compute the FirstOccurrence array and CountSymbol function
def FirstOccurrence_CountSymbol(burrows_wheeler, alphabet = ['$','A', 'C', 'G', 'T']):
    l = len(burrows_wheeler)
    CountSymbol = {}
    first_occurances = {}
    for char in alphabet:
        CountSymbol[char] = [0] * (l + 1)
    for i in range(l):
        currChar = burrows_wheeler[i]
        for char, count in CountSymbol.items():
            CountSymbol[char][i+1] = CountSymbol[char][i]
        CountSymbol[currChar][i+1] = CountSymbol[currChar][i+1] + 1
    currIndex = 0
    for char in sorted(alphabet):
        first_occurances[char] = currIndex
        currIndex = currIndex + CountSymbol[char][l]
    return first_occurances, CountSymbol


#perform pattern matching in a Burrows-Wheeler transformed string
#using the better Burrows-Wheeler matching algorithm
def BetterBWMatching(suff_arr, pattern, burrows_wheeler, starts, counts):
    occs = set()
    top = 0
    bottom = len(burrows_wheeler) - 1
    currIndex = len(pattern) - 1
    while top <= bottom:
        if currIndex >= 0:
            symbol = pattern[currIndex]
            currIndex = currIndex - 1
            if counts[symbol][bottom+1] - counts[symbol][top] > 0:
                top = starts[symbol] + counts[symbol][top]
                bottom = starts[symbol] + counts[symbol][bottom+1] - 1
            else:
                break
        else:
            for i in range(top, bottom + 1):
                occs.add(suff_arr[i])
            break

    return occs


#function to align a single read to a genome
def align_read_to_genome(read, reference_genome, k, suffix_array, burrows_wheeler, starts, counts):

    # Step 1: Break the read into k-mers
    kmers = [read[i:i+k] for i in range(len(read)-k+1)]
    #first_kmer = kmers[0]


    # Step 2: Search for matches in the BurrowsWheeler  index and extend the alignment
    best_match = None
    best_score = float('inf')
    for i, kmer in enumerate(kmers):
        positions = BetterBWMatching(suffix_array, kmer, burrows_wheeler, starts, counts)

        for pos in positions:
            # extend the alignment
            offset = i
            alignment_start = pos - offset
            alignment_end = alignment_start + len(read)
            if alignment_start < 0 or alignment_end > len(reference_genome):
                continue # alignment out of bounds, skip to next position
            ref_sequence = reference_genome[alignment_start:alignment_end]
            score = HammingDistance(read, ref_sequence)

            # check if this is the best match so far
            if score < best_score:
                best_score = score
                best_match = alignment_start

    return best_match, best_score

#function to align all reads to the genome
def align_all_reads_to_genome(suffix_array, donor_reads, reference_genome, k, burrows_wheeler, starts, counts):
    results = []
    for read_id, read_seq in donor_reads:
        best_match, best_score = align_read_to_genome(read_seq, reference_genome, k, suffix_array, burrows_wheeler, starts, counts)
        results.append({'donor_read_id': read_id,'sequence' : read_seq ,'best_match': best_match, 'best_score': best_score})
    return results


#run functions

#filename = "/content/project3a_10000_spectrum.fasta"  # Replace with the actual filename
reference_genome = GenomeAssembly(spectrum_file)
#print(reference_genome)

#Align the spectrum to sort reads in the order that they appear in the genome

k = 20

# Append the "$" symbol to the end of the genome
reconstructed_genome = reference_genome.reconstructed_genome
reconstructed_genome += "$"
#print(reconstructed_genome)

#bwt = BWT(reconstructed_genome).bwt
bwt = BurrowsWheeler(reconstructed_genome).burrows_wheeler
starts, counts = FirstOccurrence_CountSymbol(bwt)
suffix_array = suffix_array(reconstructed_genome).get_suffix_array()
donor_reads = read_fasta_file(spectrum_file)

results = align_all_reads_to_genome(suffix_array, donor_reads, reconstructed_genome, k, bwt, starts, counts)
#print(results)

# Sort the results based on the 'best_match' value
sorted_results = sorted(results, key=lambda x: x['best_match'])

'''
# Print the read IDs in the sorted order
for result in sorted_results:
    print(result['donor_read_id'])
    '''

# Open a text file for writing
with open('predictions.txt', 'w') as file:
    for result in sorted_results:
        file.write(f'>{result["donor_read_id"]}\n')
