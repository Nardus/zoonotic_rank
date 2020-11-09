# Notes:
# $(@D) means the directory of the target file


#? Usage: 'make all' updates all files as needed. To re-run the entire analysis, use 'make clean all'.
#? To change configuration options, edit the 'options.config' file in this directory
include options.config
$(info Running on $(N_CORES) cores with random seed $(RANDOM_SEED). Edit options.config to change.)
$(info )


.PHONY: all
all: download \
     get_sequences \
	 get_transcripts \
	 merge_zoonotic_status \
	 merge_and_clean_data \
	 calculate_genomic \
	 select_holdout \
	 feature_selection_runs \
	 train \
	 bag_predictions \
	 train_tax \
	 predict_novel \
	 predict_sarbeco \
	 make_plots


#? 
#? To explicitly set the path to dependencies, use make update_path
#?   This can be used multiple times
# This is needed for blast
.PHONY: update_path
update_path:
	@read -p "enter path:" path; \
	  echo "PATH=$$path:$$PATH" > .Renviron; \
		export PATH="$$path::\$$PATH"


#?
#? To run individual steps in the pipeline, combine 'make' with the command given in brackets below:

# ----------------------------------------------------------------------------------------
#?	 1. Download external data (download)
# ----------------------------------------------------------------------------------------
EXTERNALDATAFILES = ExternalData/ICTV_MasterSpeciesList_2016v1.3.xlsx \
					ExternalData/ICTV_MasterSpeciesList_2018b.xlsx \
					ExternalData/WoolhouseBrierley_2018.xlsx \
					ExternalData/Olival2017viruses.csv \
					ExternalData/Olival2017associations.csv \
					ExternalData/HousekeepingGenes.txt

.PHONY: download
download: $(EXTERNALDATAFILES)


ExternalData/ICTV_MasterSpeciesList_2016v1.3.xlsx:
	curl -L -o $@ 'https://talk.ictvonline.org/files/master-species-lists/m/msl/6776/download'
	
ExternalData/ICTV_MasterSpeciesList_2018b.xlsx:
	curl -L -o $@ 'https://talk.ictvonline.org/files/master-species-lists/m/msl/8266/download'


ExternalData/WoolhouseBrierley_2018.xlsx:
	curl -L -o $(@D)/WB2018.zip 'http://datashare.is.ed.ac.uk/download/DS_10283_2970.zip'
	unzip -u -d $(@D) $(@D)/WB2018.zip 'Woolhouse and Brierley RNA virus database.xlsx'
	mv $(@D)/'Woolhouse and Brierley RNA virus database.xlsx' $(@D)/WoolhouseBrierley_2018.xlsx
	touch $(@D)/WoolhouseBrierley_2018.xlsx   # Simply updates 'last modified' date, since unzip doesn't do this
	rm $(@D)/WB2018.zip

ExternalData/Olival2017.zip:
	curl -L -o $@ 'https://zenodo.org/record/807517/files/ecohealthalliance/HP3-v1.0.9.zip'

ExternalData/Olival2017viruses.csv: ExternalData/Olival2017.zip
	unzip -uj -d $(@D) $(@D)/Olival2017.zip 'ecohealthalliance-HP3-928327a/data/viruses.csv'
	mv $(@D)/viruses.csv $(@D)/Olival2017viruses.csv
	touch $(@D)/Olival2017viruses.csv

ExternalData/Olival2017associations.csv: ExternalData/Olival2017.zip
	unzip -uj -d $(@D) $(@D)/Olival2017.zip 'ecohealthalliance-HP3-928327a/data/associations.csv'
	mv $(@D)/associations.csv $(@D)/Olival2017associations.csv
	touch $(@D)/Olival2017associations.csv


# Gene sets:
ExternalData/HousekeepingGenes.txt:
	curl -L -o $@ 'https://www.tau.ac.il/~elieis/HKG/HK_genes.txt'



# ----------------------------------------------------------------------------------------
#?	 2. Download virus sequences from GenBank (get_sequences)
# ----------------------------------------------------------------------------------------
# Some rules below actually need the individual gb files in ExternalData/Sequences/, 
# but CombinedSequences.fasta is always the last file to be created during download,
# so if it is up to date, all sequences must be present:
ExternalData/Sequences/CombinedSequences.fasta: InternalData/Final_Accessions_Unique_Spp.csv
	python3 Misc/DownloadSequences.py

.PHONY: get_sequences
get_sequences: ExternalData/Sequences/CombinedSequences.fasta



# ----------------------------------------------------------------------------------------
#?	 3. Download human transcript sequences from Ensembl (get_transcripts)
# ----------------------------------------------------------------------------------------
CalculatedData/HumanGeneSets/TranscriptSequences.fasta: InternalData/Shaw2017_raw/ISG_CountsPerMillion_Human.csv \
                                                         ExternalData/HousekeepingGenes.txt
	python3 Misc/DownloadGeneSets.py

.PHONY: get_transcripts
get_transcripts: CalculatedData/HumanGeneSets/TranscriptSequences.fasta



# ----------------------------------------------------------------------------------------
#?	 4. Merge zoonotic status data (merge_zoonotic_status)
# ----------------------------------------------------------------------------------------
CalculatedData/ZoonoticStatus_Merged.rds: $(EXTERNALDATAFILES) \
										  InternalData/Taxonomy_UnclassifiedViruses.csv \
										  InternalData/SourcesOfZoonoses_BabayanZoonotic.csv \
										  InternalData/NameMatches_All.csv
	Rscript Scripts/MergeZoonoticStatusData.R

.PHONY: merge_zoonotic_status
merge_zoonotic_status: CalculatedData/ZoonoticStatus_Merged.rds



# ----------------------------------------------------------------------------------------
#?	 5. Merge and clean final dataset (merge_and_clean_data)
# ----------------------------------------------------------------------------------------
# This has multiple outputs: using a pattern rule ensures the command is
# run just once (see https://www.cmcrossroads.com/article/rules-multiple-outputs-gnu-make)
CalculatedData/FinalData_%.rds CalculatedData/FinalData_%.csv: InternalData/AllInternalData_Checked.csv \
															   InternalData/NameMatches_All.csv \
															   InternalData/Final_Accessions_Unique_Spp.csv \
															   CalculatedData/ZoonoticStatus_Merged.rds
	Rscript Scripts/MergeAndCleanData.R

.PHONY: merge_and_clean_data
merge_and_clean_data:  CalculatedData/FinalData_Cleaned.rds



# ----------------------------------------------------------------------------------------
#?	 6. Calculate genomic features (calculate_genomic)
# ----------------------------------------------------------------------------------------
# This actually creates multiple output files, but as they are always required
# together, simply ensuring the first of them gets updated
CalculatedData/GenomicFeatures-Virus.rds: CalculatedData/FinalData_Cleaned.rds \
                                          ExternalData/Sequences/CombinedSequences.fasta \
										  InternalData/Shaw2017_raw/ISG_PublishedData_Web.csv \
										  InternalData/Shaw2017_raw/ISG_CountsPerMillion_Human.csv \
										  CalculatedData/HumanGeneSets/TranscriptSequences.fasta
	Rscript Scripts/CalculateGenomicFeatures.R


.PHONY: calculate_genomic
calculate_genomic: CalculatedData/GenomicFeatures-Virus.rds


# ----------------------------------------------------------------------------------------
#?	 7. Remove viruses with only partial genomes available (select_holdout)
# ----------------------------------------------------------------------------------------
# NOTE: Holdout not currently used (not enough data), but this script also
#		separates out viruses with partial genomes, which should not be used 
#		for training, so we still need it
CalculatedData/%_Holdout.rds CalculatedData/%_Training.rds: CalculatedData/FinalData_Cleaned.rds
	Rscript Scripts/SelectHoldoutData.R $(RANDOM_SEED) --holdoutProportion 0


.PHONY: select_holdout
select_holdout: CalculatedData/SplitData_Holdout.rds


# ----------------------------------------------------------------------------------------
#?	 8. Compare performance of differing numbers of features (feature_selection_runs)
# ----------------------------------------------------------------------------------------
# Feature sets /  data required to calculate them:
TRAIN_REQUIREMENTS = CalculatedData/SplitData_Training.rds

TAXONOMY_REQUIREMENTS = InternalData/Taxonomy_UnclassifiedViruses.csv \
						ExternalData/ICTV_MasterSpeciesList_2018b.xlsx \

PN_REQUIREMENTS = ExternalData/Sequences/CombinedSequences.fasta

DIRECT_GENOMIC = CalculatedData/GenomicFeatures-Virus.rds
RELATIVE_GENOMIC = CalculatedData/GenomicFeatures-Distances.rds

# Data common to all possible train calls:
TRAIN_REQUIREMENTS = CalculatedData/SplitData_Training.rds


RunData/FeatureSelection_Top%: $(TRAIN_REQUIREMENTS) $(TAXONOMY_REQUIREMENTS) $(PN_REQUIREMENTS) $(DIRECT_GENOMIC) $(RELATIVE_GENOMIC)
	Rscript Scripts/TrainAndValidate.R $(RANDOM_SEED) $(notdir $(@)) --nthread $(N_CORES) \
	  --includeTaxonomy --includePN --includeVirusFeatures \
	  --includeISG --includeHousekeeping --includeRemaining \
	  --topFeatures $*

RunData/FeatureSelection_NoSelection: $(TRAIN_REQUIREMENTS) $(TAXONOMY_REQUIREMENTS) $(PN_REQUIREMENTS) $(DIRECT_GENOMIC) $(RELATIVE_GENOMIC)
	Rscript Scripts/TrainAndValidate.R $(RANDOM_SEED) $(notdir $(@)) --nthread $(N_CORES) \
	  --includeTaxonomy --includePN --includeVirusFeatures \
	  --includeISG --includeHousekeeping --includeRemaining \
	  --topFeatures 5000


SELECTION_RUN_IDS = RunData/FeatureSelection_Top10 \
					RunData/FeatureSelection_Top50 \
					RunData/FeatureSelection_Top100 \
					RunData/FeatureSelection_Top125 \
					RunData/FeatureSelection_Top150 \
					RunData/FeatureSelection_Top175 \
					RunData/FeatureSelection_Top200 \
					RunData/FeatureSelection_NoSelection

.PHONY: feature_selection_runs
feature_selection_runs: $(SELECTION_RUN_IDS)


# ----------------------------------------------------------------------------------------
#?	 9. Train models (train)
# ----------------------------------------------------------------------------------------
# NOTE: $(notdir $(@)) means the last part of the target, i.e. 'RunID' in 'RunData/RunID'

N_FEATS = 125  # Best set of models from previous step included top 125 features

#  - All possible features:
#    (already have this from previous step - FeatureSelection_Top125)


#  - Virus direct
RunData/VirusDirect: $(TRAIN_REQUIREMENTS) $(DIRECT_GENOMIC)
	Rscript Scripts/TrainAndValidate.R $(RANDOM_SEED) $(notdir $(@)) --nthread $(N_CORES) --includeVirusFeatures --topFeatures $(N_FEATS)


#  - Taxonomy
RunData/Taxonomy: $(TRAIN_REQUIREMENTS) $(TAXONOMY_REQUIREMENTS)
	Rscript Scripts/TrainAndValidate.R $(RANDOM_SEED) $(notdir $(@)) --nthread $(N_CORES) --includeTaxonomy  --topFeatures $(N_FEATS)


# - Phylogenetic neighbourhood 
#   (now using a reduced ["minimal"] set of PN features only, since we know 
#    from other work that reservoirs, etc. are not predictive)
RunData/PN: $(TRAIN_REQUIREMENTS) $(PN_REQUIREMENTS)
	Rscript Scripts/TrainAndValidate.R $(RANDOM_SEED) $(notdir $(@)) --nthread $(N_CORES) --includePN  --topFeatures $(N_FEATS)
	
RunData/PN_LongRun: $(TRAIN_REQUIREMENTS) $(PN_REQUIREMENTS)
	Rscript Scripts/TrainAndValidate.R $(RANDOM_SEED) $(notdir $(@)) --nthread $(N_CORES) --includePN  --topFeatures $(N_FEATS) \
	  --nboot 1000 --nseeds 100


#  - ISG
RunData/ISG: $(TRAIN_REQUIREMENTS) $(RELATIVE_GENOMIC)
	Rscript Scripts/TrainAndValidate.R $(RANDOM_SEED) $(notdir $(@)) --nthread $(N_CORES) --includeISG --topFeatures $(N_FEATS)

#  - Housekeeping
RunData/Housekeeping: $(TRAIN_REQUIREMENTS) $(RELATIVE_GENOMIC)
	Rscript Scripts/TrainAndValidate.R $(RANDOM_SEED) $(notdir $(@)) --nthread $(N_CORES) --includeHousekeeping --topFeatures $(N_FEATS)

#  - Remaining
RunData/Remaining: $(TRAIN_REQUIREMENTS) $(RELATIVE_GENOMIC)
	Rscript Scripts/TrainAndValidate.R $(RANDOM_SEED) $(notdir $(@)) --nthread $(N_CORES) --includeRemaining --topFeatures $(N_FEATS)


#  - All genome features
RunData/AllGenomeFeatures: $(TRAIN_REQUIREMENTS) $(RELATIVE_GENOMIC)
	Rscript Scripts/TrainAndValidate.R $(RANDOM_SEED) $(notdir $(@)) --nthread $(N_CORES) \
		--includeVirusFeatures --includeISG --includeHousekeeping --includeRemaining \
		--topFeatures $(N_FEATS)
		
RunData/AllGenomeFeatures_LongRun: $(TRAIN_REQUIREMENTS) $(RELATIVE_GENOMIC)
	Rscript Scripts/TrainAndValidate.R $(RANDOM_SEED) $(notdir $(@)) --nthread $(N_CORES) \
		--includeVirusFeatures --includeISG --includeHousekeeping --includeRemaining \
		--topFeatures $(N_FEATS) --nboot 1000 --nseeds 100



ALL_RUN_IDS = VirusDirect Taxonomy PN ISG Housekeeping Remaining \
			  AllGenomeFeatures AllGenomeFeatures_LongRun \
			  VirusDirect_ISG VirusDirect_ISG_Housekeeping PN_LongRun

TRAIN_OUTPUT_FOLDERS = $(patsubst %, RunData/%, $(ALL_RUN_IDS))

.PHONY: train
train: $(TRAIN_OUTPUT_FOLDERS)


# From here on: checking all rundata directories for the files needed / created below:
VPATH = $(TRAIN_OUTPUT_FOLDERS)


# ----------------------------------------------------------------------------------------
#?	10. Bagged predictions (bag_predictions)
# ----------------------------------------------------------------------------------------
# For long runs, increase N-models to still use top 10% of trained models:
%_LongRun_Bagged_predictions.rds: | RunData/%_LongRun
	Rscript Scripts/CalculateBaggedPredictions.R $(RANDOM_SEED) $(addsuffix _LongRun, $*) --Ntop 100

%_Bagged_predictions.rds: | RunData/%
	Rscript Scripts/CalculateBaggedPredictions.R $(RANDOM_SEED) $* --Ntop 10
 

# Currently only using bagging for long runs - need each virus to occur enough test sets:
LONG_RUN_IDS = 	AllGenomeFeatures_LongRun PN_LongRun

.PHONY: bag_predictions
bag_predictions: $(patsubst %, %_Bagged_predictions.rds, $(LONG_RUN_IDS))


# ----------------------------------------------------------------------------------------
#?	11. Fit taxonomy-based heuristic (train_tax)
# ----------------------------------------------------------------------------------------
# Get accuracy of generalising zoonotic status from family:
RunData/TaxonomyHeuristic/Test_BootstrapPredictions.rds: CalculatedData/SplitData_Training.rds
	Rscript Scripts/TrainFamilyHeuristic.R $(RANDOM_SEED) TaxonomyHeuristic --nthread $(N_CORES)


.PHONY: train_tax
train_tax: RunData/TaxonomyHeuristic/Test_BootstrapPredictions.rds


# ----------------------------------------------------------------------------------------
#?	12. Find and predict novel viruses (predict_novel)
# ----------------------------------------------------------------------------------------
# Novel viruses defined as spp added to the latest ICTV taxonomy release
#  - Matching accession numbers taken from ICTV's virus metadata resource

# Get ICTV data
ExternalData/NovelViruses/ICTV_MasterSpeciesList_2019.v1.xlsx:
	mkdir -p $(@D)
	curl -L -o $@ 'https://talk.ictvonline.org/files/master-species-lists/m/msl/9601/download'
	
ExternalData/NovelViruses/ICTV_VMR_2019.v1.xlsx:
	mkdir -p $(@D)
	curl -L -o $@ 'https://talk.ictvonline.org/taxonomy/vmr/m/vmr-file-repository/9603/download'


# Download matching sequences ($^ means all prerequisites)
ExternalData/NovelViruses/NovelViruses.gb: ExternalData/NovelViruses/ICTV_MasterSpeciesList_2019.v1.xlsx \
										   ExternalData/NovelViruses/ICTV_VMR_2019.v1.xlsx
	python3 Scripts/FindNovelViruses.py $^ ExternalData/NovelViruses/NovelViruses


# Predict
Predictions/NovelViruses.predictions.csv: ExternalData/NovelViruses/NovelViruses.gb \
										  RunData/AllGenomeFeatures_LongRun/AllGenomeFeatures_LongRun_ModelFits.rds \
										  CalculatedData/GenomicFeatures-Virus.rds \
										  CalculatedData/SplitData_Training.rds \
										  CalculatedData/GenomicFeatures-HumanCombined.rds
	Rscript Scripts/PredictNovel.R genbank ExternalData/NovelViruses/NovelViruses.gb \
								   ExternalData/NovelViruses/NovelViruses.csv \
								   Predictions/NovelViruses \
								   --random_seed $(RANDOM_SEED)


.PHONY: predict_novel
predict_novel: Predictions/NovelViruses.predictions.csv



# ----------------------------------------------------------------------------------------
#?	13. Predict Sarbecoviruses (predict_sarbeco)
# ----------------------------------------------------------------------------------------
# Get data
ExternalData/sarbecovirus/boni_et_al_NRR1_alignment.fas:
	mkdir -p $(@D)
	curl -L -o $@ 'https://raw.githubusercontent.com/plemey/SARSCoV2origins/master/alignments/sarbecovirus/NRR1/NRR1.fas'

ExternalData/sarbecovirus/sarbecovirus_raw.gb: ExternalData/sarbecovirus/boni_et_al_NRR1_alignment.fas
	python3 Scripts/get_sarbecovirus_seqs.py


# Predict
Predictions/sarbecovirus.predictions.csv: ExternalData/sarbecovirus/sarbecovirus_raw.gb \
										  RunData/AllGenomeFeatures_LongRun/AllGenomeFeatures_LongRun_ModelFits.rds \
										  CalculatedData/GenomicFeatures-Virus.rds \
										  CalculatedData/SplitData_Training.rds \
										  CalculatedData/GenomicFeatures-HumanCombined.rds
	Rscript Scripts/PredictNovel.R genbank ExternalData/sarbecovirus/sarbecovirus_raw.gb \
								   ExternalData/sarbecovirus/sarbecovirus_metadata.csv \
								   Predictions/sarbecovirus \
								   --exclude "Severe acute respiratory syndrome-related coronavirus" \
								   --random_seed $(RANDOM_SEED)


# Make a matching phylogeny for plotting
#  - This requires iqtree. Install using  conda install -c bioconda iqtree=2.0.3
#  - The phylogeny is included, so this dependency is not generally needed
#  - Using a GTR+G model, as used by Boni et al. for their ML trees (https://doi.org/10.1038/s41564-020-0771-4)
ExternalData/sarbecovirus/sarbeco_ml_phylogeny.treefile: ExternalData/sarbecovirus/boni_et_al_NRR1_alignment.fas
	iqtree -redo -s $< -m GTR+G -nt AUTO --threads-max $(N_CORES) --prefix ExternalData/sarbecovirus/sarbeco_ml_phylogeny -o "BtKY72|Bat-R_spp|Kenya|KY352407|2007-10"


.PHONY: predict_sarbeco
predict_sarbeco: Predictions/sarbecovirus.predictions.csv


# ----------------------------------------------------------------------------------------
#?	14. Plot (make_plots)
# ----------------------------------------------------------------------------------------

Plots/Figure1.pdf: CalculatedData/SplitData_Training.rds \
				   ExternalData/ICTV_MasterSpeciesList_2018b.xlsx \
				   InternalData/Taxonomy_UnclassifiedViruses.csv \
				   RunData/TaxonomyHeuristic/Test_BootstrapPredictions.rds \
				   RunData/TaxonomyHeuristic/Test_BaggedPredictions.rds \
				   $(TRAIN_OUTPUT_FOLDERS) \
				   RunData/AllGenomeFeatures_LongRun/AllGenomeFeatures_LongRun_Predictions.rds \
				   RunData/AllGenomeFeatures_LongRun/AllGenomeFeatures_LongRun_Bagged_predictions.rds
	Rscript Scripts/Plotting/MakeFigure1.R


Plots/Figure2.pdf: Plots/Figure1.pdf \
				   ExternalData/ICTV_MasterSpeciesList_2018b.xlsx \
				   CalculatedData/SplitData_Training.rds \
				   CalculatedData/GenomicFeatures-Virus.rds \
				   CalculatedData/GenomicFeatures-Distances.rds \
				   RunData/AllGenomeFeatures_LongRun/AllGenomeFeatures_LongRun_Bagged_predictions.rds
	Rscript Scripts/Plotting/MakeFigure2.R


Plots/Figure3.pdf: RunData/AllGenomeFeatures_LongRun/AllGenomeFeatures_LongRun_Bagged_predictions.rds \
				   Predictions/NovelViruses.predictions.csv \
				   Predictions/sarbecovirus.predictions.csv \
				   CalculatedData/SplitData_Training.rds \
				   Plots/Figure1.pdf \
				   ExternalData/NovelViruses/ICTV_MasterSpeciesList_2019.v1.xlsx
	Rscript Scripts/Plotting/MakeFigure3.R



## SI figures extending figure 1
# S1
Plots/Supplement_RawData.pdf: Plots/Figure1.pdf \
							  CalculatedData/SplitData_Training.rds
	Rscript Scripts/Plotting/MakeSupplementaryFigure_RawData.R


# S2
Plots/Supplement_family_auc.pdf: Plots/Figure1.pdf \
                                 CalculatedData/SplitData_Training.rds \
                                 RunData/AllGenomeFeatures_LongRun/AllGenomeFeatures_LongRun_Bagged_predictions.rds
	Rscript Scripts/Plotting/MakeSupplement_FamilyAUC.R


# S3:
Plots/Supplement_TrainingSetRanks.pdf: Plots/Figure3.pdf \
									   Plots/Figure1.pdf \
									   CalculatedData/ZoonoticStatus_Merged.rds \
									   CalculatedData/SplitData_Training.rds
	Rscript Scripts/Plotting/MakeSupplement_TrainingSetRanks.R


# S4:
Supplement_HighPriority_MissingZoonoses.pdf: Plots/Figure3.pdf \
									   		 Plots/Figure1.pdf \
											 CalculatedData/ZoonoticStatus_Merged.rds \
											 CalculatedData/SplitData_Training.rds 
	Rscript Scripts/Plotting/MakeSupplement_HighPriority_MissedZoonoses.R



## SI figures extending figure 2
# S5 & S6:
Plots/Supplement_bk_plots.pdf: Plots/Figure2.pdf \
							   ExternalData/ICTV_MasterSpeciesList_2018b.xlsx \
							   CalculatedData/SplitData_Training.rds \
							   CalculatedData/GenomicFeatures-Virus.rds
	Rscript Scripts/Plotting/MakeSupplementaryFigure_ClustersVsTaxonomy.R

Plots/Combine_tanglegrams.pdf: Plots/Supplement_bk_plots.pdf
	cd Plots && pdflatex -synctex=1 -interaction=nonstopmode Combine_tanglegrams.tex


# S7:
Plots/SupplementaryFigure_FeatureClusters.pdf: Plots/Figure2.pdf
	Rscript Scripts/Plotting/MakeSupplementaryFigure_FeatureClusters.R


# S8
Plots/SupplementaryFigure_EffectDirection.pdf: Plots/Figure2.pdf
	Rscript Scripts/Plotting/MakeSupplementaryFigure_EffectDirection.R



## SI figures extending figure 3
# S9:
Plots/Supplement_Sarbecovirus_ranks.pdf: Predictions/sarbecovirus.predictions.csv \
										 ExternalData/sarbecovirus/sarbeco_ml_phylogeny.treefile \
										 CalculatedData/SplitData_Training.rds
	Rscript Scripts/Plotting/MakeSupplement_Sarbecoviruses.R


# S10:
Plots/Supplement_methods_derived_genome_features.pdf: CalculatedData/SplitData_Training.rds \
													  Plots/Figure2.pdf \
													  CalculatedData/GenomicFeatures-Virus.rds
	Rscript Scripts/Plotting/Supplement_IllustrateDerivedGenomeFeatureCalcs.R


# S11:
Plots/Supplement_FeatureSelection.pdf: $(SELECTION_RUN_IDS)
	Rscript Scripts/Plotting/MakeSupplement_FeatureSelection.R


## SI tables
# Table SI
Plots/TableS1.csv: Plots/Figure3.pdf
	cp Plots/Intermediates/combined_virus_ranks.csv Plots/TableS1.csv


.PHONY: make_plots
make_plots: Plots/Figure1.pdf \
			Plots/Figure2.pdf \
			Plots/Figure3.pdf \
			Plots/Supplement_RawData.pdf \
			Plots/Supplement_family_auc.pdf \
			Plots/Supplement_TrainingSetRanks.pdf \
			Supplement_HighPriority_MissingZoonoses.pdf \
			Plots/Supplement_bk_plots.pdf \
			Plots/Combine_tanglegrams.pdf \
			Plots/SupplementaryFigure_FeatureClusters.pdf \
			Plots/SupplementaryFigure_EffectDirection.pdf \
			Plots/Supplement_Sarbecovirus_ranks.pdf \
			Plots/Supplement_methods_derived_genome_features.pdf \
			Plots/Supplement_FeatureSelection.pdf \
			Plots/TableS1.csv



# ----------------------------------------------------------------------------------------
# Cleanup
# ----------------------------------------------------------------------------------------
.PHONY: confirm as_distributed clean

confirm:
	@echo -n "Removing generated files - are you sure? [y/N] " && read ans && [ $${ans:-N} = y ]

#?
#? Other commands:
#?	as_distributed: Return directory to the state in which it was distributed
as_distributed: confirm
	# TODO: add a confirmation step before deleting all generated content!
	-rm -r ExternalData
	-rm -r CalculatedData
	-rm -r Plots
	-rm -r Predictions
	-rm .Renviron
	find RunData -maxdepth 1 -not -name "AllGenomeFeatures_LongRun" -not -name "PN_LongRun" -delete
	-rm RunData/AllGenomeFeatures_LongRun/AllGenomeFeatures_LongRun_Bagged_predictions.rds
	-rm RunData/AllGenomeFeatures_LongRun/AllGenomeFeatures_LongRun_Bagging_AUCs.rds
	-rm RunData/AllGenomeFeatures_LongRun/AllGenomeFeatures_LongRun_CalculatedData.rds
	-rm RunData/PN_LongRun/PN_LongRun_Bagged_predictions.rds
	-rm RunData/PN_LongRun/PN_LongRun_Bagging_AUCs.rds
	-rm RunData/PN_LongRun/PN_LongRun_CalculatedData.rds


#?	clean: Remove all intermediate files, including those required for predictions (which are distributed)
clean: as_distributed
	rm -r RunData


# ----------------------------------------------------------------------------------------
# Auto document this file (https://swcarpentry.github.io/make-novice/08-self-doc/)
#  - Comments above that start with a ? become help strings
# ----------------------------------------------------------------------------------------
.PHONY: help
help: Makefile
	@sed -n 's/^#\?//p' $<


# ----------------------------------------------------------------------------------------
# Make options
# ----------------------------------------------------------------------------------------
.DELETE_ON_ERROR:
.SECONDARY:
.NOTPARALLEL:  # Individual scripts are already parallel

#?
#?
