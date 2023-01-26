DATA_DIR = "pegLIT_kkw.csv"


def callback(data: dict):
    # tentative attribute design
    linker = peglit.pegLIT(
        seq_spacer=data["spacer"],
        seq_scaffold=data["scaffold"],
        seq_template=data["template"],
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

    batch_size = mp.cpu_count()

    # Read CSV and chunk it
    df_gen = pd.read_csv(DATA_DIR, chunksize=batch_size)
    # pooling peglit processes
    results = []
    for df_subset in tqdm(df_gen):
        data = list(df_subset.to_dict("index").values())
        with ProcessPoolExecutor(max_workers=batch_size) as executor:
            rvals = list(executor.map(callback, data))
            results.extend(
                rvals)

        # print("Done with a batch")

    pd.DataFrame.from_dict(results).rename(columns={"linker_pattern":"linker"}).to_csv(
        f"{DATA_DIR}_pegLIT_results.csv",
        index=False,
    )
