import re
import argparse
import pandas as pd
import gwaslab as gl
from pathlib import Path

class CommonVariant:
    """A class to find common variants across multiple GWAS result files and plot the results."""

    def __init__(self, input_files, pval, output, mod='_mod', weird_id=False):
        """Class constructor

        Args:
            input_files (list of str): List of input GWAS files.
            pval (float): P-value threshold for identifying common variants.
            output (str): Output file name for common variants.
            mod (str, optional): Suffix for modified SNP name files. Defaults to '_mod'.
            weird_id (bool, optional): Whether the GAWS files contain unconventional (numerical) SNP names.
                                       Default to False.
        """
        self.input_files = input_files
        self.pval = pval
        self.output = output
        self.mod = mod
        self.weird_id = weird_id

    def find_variant(self):
        """Identify common variants for a particular phenotype across all input 
        GWAS files (from different ethnicity).

        Returns:
            pd.DataFrame: DataFrame containing common variants below the p-value threshold.
        """
        dfs = []
        ethnic = []
        re_pattern = r'(eth_[^.]+)\.[^.]+\.glm\.linear$'
    
        for input_file in sorted(self.input_files):
            ethnic.append(re.search(re_pattern, input_file).group(1))
            df = pd.read_csv(input_file, sep='\t', index_col=['ID'])
            dfs.append(df['P'])
        
        p_df = pd.concat(dfs, axis=1)
        output_df = p_df[(p_df < self.pval).sum(1) == p_df.shape[1]]
        output_df.columns = ethnic
        return output_df
    
    def modify_snpname(self, file):
        """Handles GWAS output files that have non-conventional SNP ID (numeric).
        Modifies the SNP names in the file to have a specific prefix

        Args:
            file (str): Path to the input file

        Returns:
            str: Path to the modified file
        """
        df = pd.read_csv(file, sep='\t')
        df['ID'] = df['ID'].apply(lambda x: f'snp{x}')
        df.to_csv(f'{Path(file).stem}{self.mod}{Path(file).suffix}', sep='\t', index=False)
        return f'{Path(file).stem}{self.mod}{Path(file).suffix}'
    
    def plot_gwas(self, common_snp=[]):
        """Plots GWAS results with annotations for common SNPs.

        Args:
            common_snp (list of str): List of SNPs to annotate on the plot. Defaults to [].
        """
        for file in self.input_files:
            if self.weird_id:
                newfile = self.modify_snpname(file)
            else:
                newfile = file
            sumstats = gl.Sumstats(newfile, fmt='plink2')
            sumstats.plot_mqq(anno=True,
                              anno_set=common_snp,
                              pinpoint=common_snp, 
                              sig_level=self.pval,
                              stratified=True,
                              use_rank=True,
                              marker_size=(2,3),
                              save=f"{file}.png")

    def main(self):
        """Main function to find common variants and plot GWAS results.
        """
        df = self.find_variant()
        if len(df) == 0:
            return
        df.to_csv(self.output)
        common_snp = [f'snp{i}' if self.weird_id else i for i in df.index]
        self.plot_gwas(common_snp=common_snp)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find common variants and plot GWAS results.")
    parser.add_argument('-i', '--input_files', nargs='+', required=True, help='List of input GWAS files')
    parser.add_argument('-p', '--pval', type=float, required=True, help='P-value threshold')
    parser.add_argument('-o', '--output', default='common_variant.csv', help='Output file name')
    parser.add_argument('-m', '--mod', default='_mod', help='Suffix for modified SNP name files')
    parser.add_argument('-w', '--weird', default=False, help='If the GWAS contains weird SNP names')
    
    args = parser.parse_args()
    
    cv = CommonVariant(args.input_files, args.pval, args.output, args.mod, args.weird)
    cv.main()