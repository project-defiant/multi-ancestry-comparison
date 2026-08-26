import pandas as pd

CAUSAL_SNP = "snp25"


def cs_snps_susieR(path):
    df = pd.read_csv(path, sep="\t")
    return set(df.loc[df["cs_id"] == 1, "snp"])


def cs_snps_sushie(path):
    df = pd.read_csv(path, sep="\t")
    return set(df.loc[df["sushie_cs_index"] != "No CS", "snp"])


def main():
    eur_cs = cs_snps_susieR("results/EUR.susieR.cs.tsv")
    afr_cs = cs_snps_susieR("results/AFR.susieR.cs.tsv")
    sushie_cs = cs_snps_sushie("results/locus1.sushie.weights.tsv")

    print(f"EUR susieR CS size: {len(eur_cs)} -> {sorted(eur_cs)}")
    print(f"AFR susieR CS size: {len(afr_cs)} -> {sorted(afr_cs)}")
    print(f"sushie CS size: {len(sushie_cs)} -> {sorted(sushie_cs)}")

    checks = {
        "EUR CS longer than sushie CS": len(eur_cs) > len(sushie_cs),
        "AFR CS longer than sushie CS": len(afr_cs) > len(sushie_cs),
        "causal SNP in EUR CS": CAUSAL_SNP in eur_cs,
        "causal SNP in AFR CS": CAUSAL_SNP in afr_cs,
        "causal SNP in sushie CS": CAUSAL_SNP in sushie_cs,
    }
    for name, ok in checks.items():
        print(f"{'PASS' if ok else 'FAIL'}: {name}")

    if not all(checks.values()):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
