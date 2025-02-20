import pandas as pd

samples = pd.read_csv(config["deps"]["samples"], dtype=str, index_col='sample_id')
samples.index = samples.index.astype(str)
captures= pd.read_csv(config["deps"]["captures"], dtype=str, index_col='capture_id')
captures.index = captures.index.astype(str)

def lookup(df,field):
    """
    Look up field in the given dataframe
    """
    def handle_wildcards(wildcards):
        """
        Internal function for processing wildcards
        """
        result=df.loc[wildcards.capture_id, field]
        return result.split(";")
    return handle_wildcards


def get_donors_in_capture(wildcards):
    """
    Return a comma-separated list of the donors in wildcards.capture_id
    """
    capture=samples.loc[samples['donor_id']==wildcards.capture_id]
    donor_list=",".join(capture.index.tolist())
    return donor_list

def get_number_of_samples_in_capture_pool(wildcards):
    """
    Return the number of samples in wildcards.capture_id
    """
    return len(samples.loc[samples['capture_id']==wildcards.capture_id])