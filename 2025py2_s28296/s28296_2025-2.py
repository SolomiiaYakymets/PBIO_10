from io import StringIO
from Bio import Entrez, SeqIO
import pandas as pd
import matplotlib.pyplot as plt


class NCBIRetriever:
    def __init__(self, email, api_key):
        Entrez.email, Entrez.api_key = email, api_key
        Entrez.tool = 'BioScriptEx10'

    def search_taxid(self, taxid):
        try:
            with Entrez.efetch(db="taxonomy", id=taxid, retmode="xml") as h:
                org = Entrez.read(h)[0]["ScientificName"]
            print(f"Organism: {org} (TaxID: {taxid})")
            with Entrez.esearch(db="nucleotide", term=f"txid{taxid}[Organism]", usehistory="y") as h:
                res = Entrez.read(h)
            self.webenv, self.query_key = res["WebEnv"], res["QueryKey"]
            return int(res["Count"])
        except Exception as e:
            print(f"Error: {e}")
            return 0

    def fetch_records(self, start=0, max_records=10):
        if not hasattr(self, 'webenv') or not hasattr(self, 'query_key'):
            return []
        try:
            with Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", retstart=start,
                               retmax=min(max_records, 500), webenv=self.webenv, query_key=self.query_key) as h:
                return h.read()
        except Exception as e:
            print(f"Error: {e}")
            return ""


def parse_and_filter(text, min_len, max_len):
    records = []
    for r in SeqIO.parse(StringIO(text), "genbank"):
        l = len(r.seq)
        if min_len <= l <= max_len:
            records.append({"Accession": r.id, "Length": l, "Description": r.description})
    return records


def save_csv(records, filename):
    pd.DataFrame(records).to_csv(filename, index=False)


def save_chart(records, filename):
    df = pd.DataFrame(records).sort_values(by="Length", ascending=False)
    plt.figure(figsize=(12, 6))
    plt.plot(df["Accession"], df["Length"], marker='o')
    plt.xticks(rotation=90, fontsize=8)
    plt.xlabel("GenBank Accession")
    plt.ylabel("Sequence Length")
    plt.title("GenBank Sequence Lengths")
    plt.tight_layout()
    plt.savefig(filename)


def main():
    email = input("Enter your email address for NCBI: ")
    api_key = input("Enter your NCBI API key: ")
    taxid = input("Enter taxid: ")
    min_len, max_len = int(input("Min length: ")), int(input("Max length: "))
    r = NCBIRetriever(email, api_key)
    if r.search_taxid(taxid) == 0:
        return
    data = r.fetch_records(max_records=100)
    filtered = parse_and_filter(data, min_len, max_len)
    if not filtered:
        return
    with open(f"taxid_{taxid}_sample.gb", "w") as f:
        f.write(data)
    save_csv(filtered, f"taxid_{taxid}_filtered.csv")
    save_chart(filtered, f"taxid_{taxid}_chart.png")
    print("All files saved.")


if __name__ == "__main__":
    main()
