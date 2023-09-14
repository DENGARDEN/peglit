DATA_DIR = "./4k_small_lib_epeg-2023-09-13_22-44-10.csv"
SPCAS9_SCAFFOLD = (
    "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC"
)

OPTIMIZED_SGRNA_SCAFFOLD = "GTTTcAGAGCTAtgctgGAAAcagcaTAGCAAGTTgAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC"  # doi://10.1186/s13059-015-0846-3 Fig. 2b T2C-4
tevopreQ1_DP = "CGCGGTTCTATCTAGTTACGCGTTAAACCAACTAGAA"
LINKER_PATTERN = "NNNNNNNN"


def callback(data: dict):
    # tentative attribute design
    linker = peglit.pegLIT(
        seq_spacer=data["Guide"],
        seq_scaffold=data["Scaffold"],
        seq_template=data["RTT"],
        seq_pbs=data["PBS"],
        linker_pattern=data["linker_pattern"],
        seq_motif=data["motif"],
    )
    data["linker_pattern"], *_ = linker
    return data


if __name__ == "__main__":
    import multiprocessing as mp
    from concurrent.futures import ProcessPoolExecutor

    import peglit
    import pandas as pd
    from tqdm import tqdm

    tqdm.pandas()
    batch_size = mp.cpu_count()
    # batch_size =1

    # Read CSV and chunk it
    # Some modification for pe6 case
    df = pd.read_csv(DATA_DIR)
    # df["linker_pattern"] = str(LINKER_PATTERN)
    # df["motif"] = tevopreQ1_DP
    df["Scaffold"] = OPTIMIZED_SGRNA_SCAFFOLD
    df.to_csv(DATA_DIR)
    del df

    df_gen = pd.read_csv(DATA_DIR, chunksize=batch_size)
    # pooling peglit processes
    results = []
    for df_subset in tqdm(df_gen):
        data = list(df_subset.to_dict("index").values())
        with ProcessPoolExecutor(max_workers=batch_size) as executor:
            rvals = list(executor.map(callback, data))
            results.extend(rvals)

        # print("Done with a batch")

    pd.DataFrame.from_dict(results).rename(columns={"linker_pattern": "linker"}).to_csv(
        f"{DATA_DIR}_pegLIT_results.csv",
        index=False,
    )
