rule aggregate_souporcell:
    """
    This rule does nothing interesting, but can be used to trigger a rerun of just the souporcell branch of the DAG
    """
    input:
        all_souporcell_results=expand(
                out_dir/"souporcell/annotation/vireosnp/{capture_id}/clusters_annotated.tsv",
                capture_id=captures.index.tolist()
        ),
        # de-novo-vs-known partition QC depends on the known-GT second mode, so
        # make its outputs reachable via the souporcell target too.
        all_known_results=expand(
                out_dir/"souporcell/{capture_id}/known/clusters.tsv",
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
                out_dir/"vireo/{capture_id}/GT_donors.vireo.vcf.gz",
                out_dir/"vireo/{capture_id}/ambient_summary.csv",
                out_dir/"vireo/{capture_id}/ambient_stats.csv",
            ], 
            capture_id=captures.index.tolist()
        )
    output:
        checkpoint_file=out_dir/"vireo/all.done"
    shell:
        """
        echo {input.all_vireo_results} > {output.checkpoint_file}
        """

# NOTE on selectable tools (config params.tools): the four per-tool sentinels
# above (aggregate_souporcell/demuxlet/freemuxlet/vireo) are only ever pulled
# into the DAG by `rule all`, which now requires each one ONLY when that tool is
# enabled -- so a disabled tool's branch is never built. aggregate_demuxafy
# depends solely on the tool-agnostic cell_assignment.tsv (produced by the
# demuxafy->wrangle chain), so it does NOT demand demuxlet/freemuxlet outputs;
# the demuxafy rule itself declares only the enabled tools' outputs as inputs.
rule aggregate_demuxafy:
    input:
        expand(out_dir/"assignment/{capture_id}/cell_assignment.tsv", capture_id=captures.index.tolist())
    output:
        checkpoint_file=out_dir/"assignment/all.done"
    shell:
        """
        echo {input} > {output.checkpoint_file}
        """

rule aggregate_hier:
    """
    Sentinel for the hierarchical (sex-aware) consensus branch, run alongside
    demuxafy. Collects the per-capture hierarchical_consensus.tsv tables.
    """
    input:
        expand(out_dir/"assignment/{capture_id}/hierarchical_consensus.tsv", capture_id=captures.index.tolist())
    output:
        checkpoint_file=out_dir/"assignment/hierarchical.done"
    shell:
        """
        echo {input} > {output.checkpoint_file}
        """