rule aggregate_souporcell:
    """
    This rule does nothing interesting, but can be used to trigger a rerun of just the souporcell branch of the DAG
    """
    input:
        all_souporcell_results=expand(
                out_dir/"souporcell/annotation/vireosnp/{capture_id}/clusters_annotated.tsv", 
                capture_id=captures.index.tolist()
        )
    output:
        checkpoint_file=out_dir/"souporcell/all.done"
    shell:
        """
        echo {input.all_souporcell_results} > {output.checkpoint_file}
        """