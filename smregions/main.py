"""
Contains the command line parsing
"""

from os import path

import click
import bglogs

from smregions import config
from smregions.smregions import SMRegions


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


def main(mutations_file, elements_file, regions_file, signature_file, output_file, config_file, config_override_dict=None):

    # Skip if done
    if path.exists(output_file):
        bglogs.warning("Already calculated at '{}'".format(output_file))
        return
    
    configuration = config.load(config_file, override=config_override_dict)

    analysis = SMRegions(mutations_file, elements_file, regions_file, signature_file, output_file, configuration)

    bglogs.info('Running analysis')
    # Run the analysis
    analysis.run()


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-m', '--muts', 'mutations_file', type=click.Path(exists=True), help='Variants file', metavar='MUTATIONS_FILE',required=True)
@click.option('-e', '--elements', 'elements_file', type=click.Path(exists=True), metavar='ELEMENTS_FILE', help='Genomic elements to analyse', required=True)
@click.option('-r', '--regions', 'regions_file', type=click.Path(exists=True), metavar='REGIONS_FILE', help='Genomic regions of interest', required=True)
@click.option('-s', '--signature', 'signature_file', type=click.Path(exists=True), metavar='SIGNATURE_FILE', help='Signature file. Default equial probabilities', default=None)
@click.option('-o', '--output', 'output_folder', type=click.Path(), metavar='OUTPUT_FOLDER', help="Output folder. Default to regions file name without extensions.", default=None)
@click.option('-c', '--configuration', 'config_file', default=None, type=click.Path(exists=True), metavar='CONFIG_FILE', help="Configuration file. Default to 'smregions.conf' in the current folder if exists or to ~/.config/bbglab/smregions.conf if not.")
@click.option('--cores', default=1)
@click.option('--debug', help="Show more progress details", is_flag=True)
@click.version_option()
def cmdline(mutations_file, elements_file, regions_file, signature_file, output_folder, config_file, cores, debug):
    """
    Run SMRegions on the genomic regions in ELEMENTS FILE and the regions of interest REGIONS_FILE
    using the mutations in MUTATIONS FILE.

    """
    bglogs.configure(debug=True if debug else False)

    override_config = {'cores': cores}

    main(mutations_file, elements_file, regions_file, signature_file, output_folder, config_file, override_config)


if __name__ == "__main__":
    cmdline()
