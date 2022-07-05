#
# Utility functions for parsing coding sequences from coordinates
#
require(seqinr)

check_cds <- function(sequence) {
	# Check a cds for validity:
	# - Length must be divisible by 3
	# - Translation should not contain stop codons
	if (length(sequence) %% 3 != 0)
		stop("Length of extracted CDS is not a multiple of 3. Check coordinates.")
	
	translation <- seqinr::getTrans(sequence)
	
	if ("*" %in% translation[1:(length(translation) - 1)])
		stop("Extracted CDS contained an internal stop codon.")
}


extract_cds_internal <- function(sequence, start_coordinate, stop_coordinate, allow_complementary) {
	reverse_complement <- FALSE
	
	# Check inputs:
	if (!is.SeqFastadna(sequence))
		stop("Invalid sequence input: 'sequence' should be a single sequence as produced by the seqinr library (a SeqFastadna object).")
	
	if (length(sequence) == 1) 
		stop("Sequence should be a vector: Ensure 'as.string' is set to FALSE when using read.fasta().")
	
	if (start_coordinate >= stop_coordinate) {
		if (allow_complementary) {
			warning("Start coordinate occurs after stop coordinate - assuming open reading frame is on the complementary strand.")
			reverse_complement <- TRUE
			tempstop <- start_coordinate
			start_coordinate <- stop_coordinate
			stop_coordinate <- tempstop
			
		} else {
			stop("Start coordinate must be lower than stop coordinate.")
		}
	}
	
	if (start_coordinate == 0)
		warning("Potential 0-based coordinate detected. The first position is at coordinate 1.")
	
	# Extract coding sequence:
	cds <- sequence[start_coordinate:stop_coordinate]
	cds <- seqinr::as.SeqFastadna(cds, name = attr(sequence, "name"), Annot = attr(sequence, "Annot"))
	
	# Reverse complement (if needed)
	if (reverse_complement) {
		cds <- rev(seqinr::comp(cds))
	}
	
	# Check cds for validity before returning
	check_cds(cds)
	
	return(cds)
}

# User-facing version: capture error messages and report the exact CDS which was problematic:
extract_cds <- function(sequence, start_coordinate, stop_coordinate, allow_complementary = FALSE) {
	tryCatch(
		extract_cds_internal(sequence, start_coordinate, stop_coordinate, allow_complementary),
		error = function(error_message) {
			outMessage <- sprintf("Extracting CDS failed (sequence '%s', position %d to %d), with message:\n%s",
								  attr(sequence, "name"), start_coordinate, stop_coordinate, error_message)
			stop(outMessage)
		}
	)
}