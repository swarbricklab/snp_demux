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

rule aggregate_demuxlet:
    """
    This rule does nothing interesting, but can be used to trigger a rerun of just the demuxlet branch of the DAG
    """
    input:
        all_demuxlet_results=expand(
                                out_dir/"demuxlet/{capture_id}/summary.txt", 
                                capture_id=captures.index.tolist()
        )
    output:
        checkpoint_file=out_dir/"demuxlet/all.done"
    shell:
        """
        echo {input.all_demuxlet_results} > {output.checkpoint_file}
        """

rule aggregate_freemuxlet:
    """
    This rule does nothing interesting, but can be used to trigger a rerun of just the freemuxlet branch of the DAG
    """
    input:
        all_freemuxlet_results=expand(
                                out_dir/"freemuxlet/{capture_id}/summary.txt", 
                                capture_id=captures.index.tolist()
        )
    output:
        checkpoint_file=out_dir/"freemuxlet/all.done"
    shell:
        """
        echo {input.all_freemuxlet_results} > {output.checkpoint_file}
        """

rule aggregate_vireo:
    """
    This rule does nothing interesting, but can be used to trigger a rerun of just the demuxlet branch of the DAG
    """
    input:
        all_vireo_results=expand(
            [
                out_dir/"vireo/{capture_id}/summary.tsv", 
                out_dir/"vireo/{capture_id}/GT_donors.vireo.vcf.gz"
            ], 
            capture_id=captures.index.tolist()
        )
    output:
        checkpoint_file=out_dir/"vireo/all.done"
    shell:
        """
        echo {input.all_vireo_results} > {output.checkpoint_file}
        """

rule aggregate_demuxafy:
    input:
        expand(out_dir/"assignment/{capture_id}/cell_assignment.tsv", capture_id=captures.index.tolist())
    output:
        checkpoint_file=out_dir/"assignment/all.done"
    shell:
        """
        echo {input} > {output.checkpoint_file}
        """