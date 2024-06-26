
# import glob
# import pandas as pd


# remove white space and fix parentethies in tissue name
fixChr = str.maketrans({" ": None, "(": "_", ")": "_"})

# Tissues with multiple sub tissues types, names are sub-tissue type names
# include Kidney-Cortex
GtexMultiSubTissues = 'Adipose-Subcutaneous,Adipose-Visceral(Omentum),WholeBlood,Cells-EBV-transformedlymphocytes,Artery-Tibial,Artery-Aorta,Artery-Coronary,Brain-Cortex,Brain-Caudate(basalganglia),Brain-Nucleusaccumbens(basalganglia),Brain-Cerebellum,Brain-CerebellarHemisphere,Brain-FrontalCortex(BA9),Brain-Putamen(basalganglia),Brain-Hypothalamus,Brain-Hippocampus,Brain-Anteriorcingulatecortex(BA24),Brain-Spinalcord(cervicalc-1),Brain-Amygdala,Brain-Substantianigra,Colon-Transverse,Colon-Sigmoid,Esophagus-Mucosa,Esophagus-Muscularis,Esophagus-GastroesophagealJunction,Heart-LeftVentricle,Heart-AtrialAppendage,Kidney-Cortex,Skin-SunExposed(Lowerleg),Skin-NotSunExposed(Suprapubic),Cells-Culturedfibroblasts'

GtexSingleSubTissues = 'Muscle-Skeletal,Thyroid,Nerve-Tibial,Lung,Breast-MammaryTissue,Testis,Stomach,Pancreas,Pituitary,AdrenalGland,Prostate,Spleen,Liver,SmallIntestine-TerminalIleum,Ovary,MinorSalivaryGland,Vagina,Uterus,Bladder'
GtexOutputSubTissues = GtexSingleSubTissues + ',' + GtexMultiSubTissues
GtexOutputSubTissues = GtexOutputSubTissues.translate(fixChr)

print(f"Number of Tissues: {len(GtexOutputSubTissues.split(','))}")
print('\n'.join(GtexOutputSubTissues.split(',')))


localrules: 
  all

rule all:
    input: 
        expand('smk-plots/gtex-sqtl-enrichment/{TISSUE}-qqplot.png', TISSUE=GtexOutputSubTissues.split(',')),


rule plotGTExsQTLEnrichment:
    output:
        qqplot = 'smk-plots/gtex-sqtl-enrichment/{tissue}-qqplot.png',
        scatter = 'smk-plots/gtex-sqtl-enrichment/{tissue}-scatter.png',
        rds = 'smk-plots/gtex-sqtl-enrichment/{tissue}-plotdata.rds',
    params:
        script = 'scripts/plot-gtex-sqtl-enrichment.R',
        basepath = '/project/yangili1/cdai/SpliFi',
        out_prefix = '/project/yangili1/cdai/splice-pub/smk-plots/gtex-sqtl-enrichment',
        FDR = 0.1
    shell:
        '''
        Rscript {params.script} {wildcards.tissue} {params.basepath} {params.out_prefix} {params.FDR}
        touch {output}
        '''

