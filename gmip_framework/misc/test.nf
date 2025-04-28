nextflow.enable.dsl=2

workflow {
    chr_list = (1..22).collect()

    channel_of_chr = Channel.from(chr_list)

    process_echo_chromosomes(channel_of_chr)
}

process process_echo_chromosomes {
    input:
    val chr

    script:
    """
    echo "Chromosome $chr"
    """
}
