
# import glob
# import pandas as pd


# remove white space and fix parentethies in tissue name
fixChr = str.maketrans({" ": None, "(": "_", ")": "_"})

# Tissues with multiple sub tissues types, names are sub-tissue type names
# include Kidney-Cortex
GtexMultiSubTissues = 'Adipose-Subcutaneous,Adipose-Visceral(Omentum),WholeBlood,Cells-EBV-transformedlymphocytes,Artery-Tibial,Artery-Aorta,Artery-Coronary,Brain-Cortex,Brain-Caudate(basalganglia),Brain-Nucleusaccumbens(basalganglia),Brain-Cerebellum,Brain-CerebellarHemisphere,Brain-FrontalCortex(BA9),Brain-Putamen(basalganglia),Brain-Hypothalamus,Brain-Hippocampus,Brain-Anteriorcingulatecortex(BA24),Brain-Spinalcord(cervicalc-1),Brain-Amygdala,Brain-Substantianigra,Colon-Transverse,Colon-Sigmoid,Esophagus-Mucosa,Esophagus-Muscularis,Esophagus-GastroesophagealJunction,Heart-LeftVentricle,Heart-AtrialAppendage,Kidney-Cortex,Skin-SunExposed(Lowerleg),Skin-NotSunExposed(Suprapubic),Cells-Culturedfibroblasts'

GtexSingleSubTissues = 'Muscle-Skeletal,Thyroid,Nerve-Tibial,Lung,Breast-MammaryTissue,Testis,Stomach,Pancreas,Pituitary,AdrenalGland,Prostate,Spleen,Liver,SmallIntestine-TerminalIleum,Ovary,MinorSalivaryGland,Vagina,Uterus'
# GtexSingleSubTissues = 'Muscle-Skeletal,Thyroid,Nerve-Tibial,Lung,Breast-MammaryTissue,Testis,Stomach,Pancreas,Pituitary,AdrenalGland,Prostate,Spleen,Liver,SmallIntestine-TerminalIleum,Ovary,MinorSalivaryGland,Vagina,Uterus,Bladder'
GtexOutputSubTissues = GtexSingleSubTissues + ',' + GtexMultiSubTissues
GtexOutputSubTissues = GtexOutputSubTissues.translate(fixChr)


print(f"Number of Tissues: {len(GtexOutputSubTissues.split(','))}")
# print('\n'.join(GtexOutputSubTissues.split(',')))


localrules: 
  all

rule all:
  input: 
    expand('smk-plots/gtex-sqtl-enrichment-v4/{TISSUE}-qqplot.svg', TISSUE=GtexOutputSubTissues.split(',')),
    # expand('smk-plots/gtex-sqtl-enrichment-v4-C/{TISSUE}-qqplot.svg', TISSUE=GtexOutputSubTissues.split(',')),
    expand('data/torus/torus_{TISSUE}_enrichment_p-sqtl.txt', TISSUE=GtexOutputSubTissues.split(',')),


rule plotGTExsQTLEnrichment:
  # previous version was saved smk-plots/gtex-sqtl-enrichment-v1
  output:
    qqplot = 'smk-plots/gtex-sqtl-enrichment-v4/{tissue}-qqplot.svg',
    scatter = 'smk-plots/gtex-sqtl-enrichment-v4/{tissue}-scatter.svg',
    rds = 'smk-plots/gtex-sqtl-enrichment-v4/{tissue}.rds',
  params:
    script = 'scripts/plot-gtex-sqtl-enrichment.R',
    basepath = '/project/yangili1/cdai/SpliFi',
    out_prefix = '/project/yangili1/cdai/splice-pub/smk-plots/gtex-sqtl-enrichment-v4',
    FDR = 0.1,
    annot = '/project/yangili1/cdai/splice-pub/data/WGS_Feature_overlap_collapsed_VEP_short_4torus.MAF01.txt.gz',
  shell:
    '''
    Rscript {params.script} {wildcards.tissue} {params.basepath} {params.out_prefix} {params.FDR} {params.annot}
    touch {output}
    '''

use rule plotGTExsQTLEnrichment as plotGTExsQTLEnrichment_test with:
  output:
    qqplot = 'smk-plots/gtex-sqtl-enrichment-v4-C/{tissue}-qqplot.svg',
    scatter = 'smk-plots/gtex-sqtl-enrichment-v4-C/{tissue}-scatter.svg',
    rds = 'smk-plots/gtex-sqtl-enrichment-v4-C/{tissue}.rds',
  params:
    script = 'scripts/plot-gtex-sqtl-enrichment-test.R',
    basepath = '/project/yangili1/cdai/SpliFi', out_prefix = '/project/yangili1/cdai/splice-pub/smk-plots/gtex-sqtl-enrichment-v4-C',
    FDR = 0.1,
    annot = '/project/yangili1/cdai/splice-pub/data/WGS_Feature_overlap_collapsed_VEP_short_4torus.MAF01.txt.gz',


rule prepTorusInputData:
  output: 
    prod = 'data/torus/torus_{tissue}_p-sqtl.txt.gz',
    unpr = 'data/torus/torus_{tissue}_u-sqtl.txt.gz',
    snpmap = 'data/torus/torus_{tissue}_snp_map.txt.gz',
    genemap = 'data/torus/torus_{tissue}_gene_map.txt.gz',
  params:
    script = 'scripts/torus-make-input-data.py',
    fdr = 0.1,
    outprefix = 'data/torus/torus'
  shell:
    '''
    python {params.script} --tissue {wildcards.tissue} --fdr {params.fdr} --outprefix {params.outprefix}
    touch {output}

    '''



# torus annotation columns (18)
# but some columns dont work
# SNP enhancer_d promoter_d open_chromatin_region_d promoter_flanking_region_d CTCF_binding_site_d TF_binding_site_d 3_prime_UTR_variant_d 5_prime_UTR_variant_d frameshift_variant_d intron_variant_d missense_variant_d non_coding_transcript_exon_variant_d splice_acceptor_variant_d splice_donor_variant_d splice_region_variant_d stop_gained_d synonymous_variant_d
rule RunTorusEnrich:
  input:
    prod = 'data/torus/torus_{tissue}_p-sqtl.txt.gz',
    unpr = 'data/torus/torus_{tissue}_u-sqtl.txt.gz',
    snpmap = 'data/torus/torus_{tissue}_snp_map.txt.gz',
    genemap = 'data/torus/torus_{tissue}_gene_map.txt.gz',
    annot = 'data/WGS_Feature_overlap_collapsed_VEP_short_4torus.MAF01.space.txt.gz'
  output:
    p_enri = 'data/torus/torus_{tissue}_enrichment_p-sqtl.txt',
    u_enri = 'data/torus/torus_{tissue}_enrichment_u-sqtl.txt',
  log: 'logs/torus_{tissue}_enrichment.log'
  shell:
    '''
    module load boost gsl
    torus -est -d {input.prod} -smap {input.snpmap} -gmap {input.genemap} -annot <(zcat {input.annot} | cut -d " " -f 1-9 | bgzip -c) 1> {output.p_enri} 2> {log}

    torus -est -d {input.unpr} -smap {input.snpmap} -gmap {input.genemap} -annot <(zcat {input.annot} | cut -d " " -f 1-9 | bgzip -c) 1> {output.u_enri} 2>> {log}
    '''

