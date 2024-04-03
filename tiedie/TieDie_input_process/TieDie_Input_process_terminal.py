import argparse
import omnipath as op
import pandas as pd
import numpy as np
import csv

def parse_arguments():
    parser = argparse.ArgumentParser(description='Process input files.')
    parser.add_argument('-tf','--transcriptomics_file', required=True, help='Path to transcriptomics file')
    parser.add_argument("--value_column", help="Column number for Z-score values", dest="value_column", action="store", type=int, required=True)
    parser.add_argument('--endpoint_file', required=True, help='Path to endpoint file')
    parser.add_argument('--hmi_prediction_file', required=True, help='Path to hmi prediction output file')
    parser.add_argument('--output_dir', required=True, help='Path to the output directory (resource)')
    parser.add_argument("--upstream_input_filename", type=str, default="upstream.input", help="Custom name for upstream output file")
    parser.add_argument("--downstream_input_filename", type=str, default="downstream.input", help="Custom name for downstream output file")
    parser.add_argument("--pathway_input_filename", type=str, default="pathway.sif", help="Custom name for pathway output file")
    
    return parser.parse_args()

def main():
    args = parse_arguments()

    ppis = op.interactions.OmniPath.get(genesymbols=1)
    
    tf_tg = op.interactions.Transcriptional.get(databases='CollecTRI', genesymbols=1)
    #print(tf_tg.columns)
    with open(args.transcriptomics_file) as transcriptomics:
        genes = []
        transcriptomics.readline()
        for line in transcriptomics:
            line = line.strip().split(',')

            #Python starts to count by 0; therefore the user-provided column number should be decreased by one
            value_column_new = int(args.value_column) - 1
            if float(line[value_column_new]):
                genes.append(line[0])

    # Creating the contextualised ppi DataFrame
    contextulised_ppi = ppis[ppis['source_genesymbol'].isin(genes) & ppis['target_genesymbol'].isin(genes)]

    # Filtering endpoint genes
    endpoint_genes = pd.read_csv(args.endpoint_file, sep=',')
    endpoint_genes = endpoint_genes[endpoint_genes['FDR'] < 0.05]
    endpoint_genes['gene_name']

    # Filtering tf_tg based on conditions
    contextulised_tf_tg = tf_tg[tf_tg['source_genesymbol'].isin(genes) & tf_tg['target_genesymbol'].isin(endpoint_genes['gene_name'])]
    contextulised_tf_tg = contextulised_tf_tg.drop_duplicates()

    # Get list of regulators with the number of targeted genes
    regs = contextulised_tf_tg.pivot_table(columns=['source'], aggfunc='size')
    regs = pd.DataFrame(regs)
    regs = regs.reset_index()
    regs = regs.rename(columns={"source": "TF_name", 0: "num_degs"})

    # Save intermediate results
    contextulised_tf_tg.to_csv(f"{args.output_dir}/contextualised_regulator-deg_network.txt", sep="\t", index=False)
    regs.to_csv(f"{args.output_dir}/contextualised_regulators_of_degs.txt", sep="\t", index=False)


    # Extract relevant columns from contextualised PPI
    ppis2 = contextulised_ppi[['source', 'consensus_stimulation', 'target']].copy()

    # Update consensus_stimulation values
    ppis2.loc[ppis2['consensus_stimulation'] == 1, 'consensus_stimulation'] = 'stimulates>'
    ppis2.loc[ppis2['consensus_stimulation'] == 0, 'consensus_stimulation'] = 'inhibits>'

    # Drop duplicates and rename columns
    ppis2 = ppis2.drop_duplicates()
    ppis2 = ppis2.rename(columns={'consensus_stimulation': 'direction'})

    # Save pathway.sif file
    ppis2.to_csv(f"{args.output_dir}/{args.pathway_input_filename}", sep='\t', index=None, header=False)

    # Load bacterial-human binding protein interactions
    hbps = pd.read_csv(args.hmi_prediction_file, sep="\t", index_col=False)

    # Extract relevant columns and set default sign value
    hbps = hbps[['# Human Protein', 'Bacterial protein']].copy()
    hbps['sign'] = '+'

    # Process upstream input
    if 'sign' in hbps.columns:
        hbps_f = hbps[(hbps['sign'] == '-') | (hbps['sign'] == '+')]
        if len(hbps) != len(hbps_f):
            print("WARNING: Some of the bacterial-human binding protein interactions were discarded as the values in the 'sign' column were not '+' or '-'.")
    
        if len(hbps_f) == 0:
            print("ERROR: bacterial-human binding protein interactions do not have the correct values in 'sign' column. They should be '+' or '-'.")
    
        hbps2 = hbps[['# Human Protein', 'sign']].copy()
        hbps2 = hbps2.rename(columns={'sign': 'direction'})
        hbps2 = hbps2.groupby(['# Human Protein', 'direction']).size().reset_index(name='n')
        hbps2 = hbps2[['# Human Protein', 'n', 'direction']]
    else:
        hbps2 = hbps[['# Human Protein']].copy()
        hbps2 = hbps2.groupby('# Human Protein').size().reset_index(name='n')
        hbps2['direction'] = '-'

    # Save upstream input file
    hbps2.to_csv(f"{args.output_dir}/{args.upstream_input_filename}", sep='\t', index=None, header=False)
    
    # Create bacterial annotation
    bac_annotation = hbps[['Bacterial protein','Bacterial protein']].copy()
    bac_annotation.reset_index(drop=True, inplace=True)

    # Process downstream input
    genes = endpoint_genes
    genes = genes.rename(columns={'gene_name':'target_genesymbol'})

    # Join the tf-deg network with the deg lfc values
    tfs = contextulised_tf_tg
    tfs2 = tfs.merge(genes, on='target_genesymbol')
    print(tfs2.columns)
    tfs2_filtered = tfs2[['source', 'target', 'consensus_stimulation', 'logFC']].copy()
    tfs2_filtered['exp_sign'] = np.where(tfs2_filtered['consensus_stimulation'] == True, tfs2_filtered['logFC'], -1*(tfs2_filtered['logFC']))

    # Get the sum of all lfc*sign values - and the number of target genes for each tf
    tfs3 = tfs2_filtered[['source', 'target', 'exp_sign']].copy()
    tfs3_1 = tfs3.groupby('source')
    tfs4 = tfs3_1.count()
    tfs4['sumof'] = tfs3_1.exp_sign.sum()
    tfs4['final_val'] = tfs4['sumof'] / tfs4['exp_sign']
    tfs4['sign'] = np.where(tfs4['final_val'] >= 0, '+', '-')
    tfs4 = tfs4.drop(['exp_sign', 'sumof', 'target'], axis=1)

    # Save downstream data with custom name
    downstream_filename = args.downstream_input_filename
    tfs4.to_csv(f"{args.output_dir}/{args.downstream_input_filename}", sep='\t', header=False)


if __name__ == "__main__":
    main()
