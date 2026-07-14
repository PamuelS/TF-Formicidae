MOTIFS = glob_wildcards("../../03_tf_research/05_aggregate/02_totalscore/totalscore_{motif}.tsv").motif
PHENOTYPES = ["pgls1_parasitism", "pgls2_larval_feeding", "pgls3_cocoon", "pgls4a_wings_ergatoid", "pgls4b_wings_brachypterous", "pgls5_castes", "pgls6a_queens_binary", "pgls6b_queens_ordered"]

rule all:
	input:
		expand("../02_pgls_Ortho_results/{motif}/pgls5_Ortho_castes.csv", motif=MOTIFS, phenotype=PHENOTYPES)

rule pgls_phenotypes:
	input:
		tree = "tree_rooted.nwk",
		phenotypes = "updated_dataset.tsv",
		motifs     = "../../03_tf_research/05_aggregate/05_totalscore_orthofinder/totalscore_{motif}.tsv"
	output:
		pgls5  = "../02_pgls_Ortho_results/{motif}/pgls5_Ortho_castes.csv"
	log:
		"../01_pgls_results/{motif}/logs/pgls_phenotypes.log"
	params:
		cores = "40"
	shell:
		# alternative way to run the script so that the parallelization can work properly setting the arguments for the Rscript
		"OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 Rscript pgls_Ortho_phenothypes.R {input.tree} {input.phenotypes} {input.motifs} {params.cores} {output.pgls5}"
