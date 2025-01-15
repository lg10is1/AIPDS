import pandas as pd
import os
from Bio.Blast.Applications import NcbiblastpCommandline

# 转换CSV文件为FASTA格式
def csv_to_fasta(csv_file, fasta_file):
    df = pd.read_csv(csv_file)
    with open(fasta_file, 'w') as f:
        for i, sequence in enumerate(df["Heavy Chain Sequence"]):
            f.write(f">Seq_{i+1}\n{sequence}\n")
    print(f"FASTA 文件已保存为 {fasta_file}")

# 转换TXT文件为FASTA格式
def txt_to_fasta(txt_file, fasta_file):
    with open(txt_file, 'r') as f, open(fasta_file, 'w') as out:
        for i, line in enumerate(f.readlines()):
            out.write(f">NewSeq_{i+1}\n{line.strip()}\n")
    print(f"新序列的 FASTA 文件已保存为 {fasta_file}")

# 运行 BLAST 并输出结果
def run_blast(query_fasta, db, output_file, evalue=0.001):
    blastp_cline = NcbiblastpCommandline(
        query=query_fasta,
        db=db,
        evalue=evalue,
        outfmt=6,  # 输出格式为表格（Tabular）
        out=output_file
    )
    stdout, stderr = blastp_cline()
    print(f"BLAST 分析完成，结果已保存为 {output_file}")

# 解析 BLAST 输出结果为相似性矩阵
def parse_blast_output(output_file, query_fasta, db_fasta):
    # 获取查询和数据库序列的 ID 列表
    query_ids = [line[1:].strip() for line in open(query_fasta) if line.startswith(">")]
    db_ids = [line[1:].strip() for line in open(db_fasta) if line.startswith(">")]

    # 初始化相似性矩阵
    similarity_matrix = pd.DataFrame(0, index=query_ids, columns=db_ids)

    # 填充矩阵
    with open(output_file, 'r') as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                query_id, subject_id, identity = parts[0], parts[1], float(parts[2])
                if query_id in similarity_matrix.index and subject_id in similarity_matrix.columns:
                    similarity_matrix.loc[query_id, subject_id] = identity

    return similarity_matrix

# 主函数
def main():
    # 文件路径
    csv_file = "sabdab_cdr_seq_label_all.csv"
    txt_file = "new_sequences.txt"
    db_fasta = "database.fasta"
    query_fasta = "new_sequences.fasta"
    output_file = "blast_results.txt"
    matrix_output_file = "similarity_matrix_blast.csv"

    # 转换 CSV 和 TXT 为 FASTA
    csv_to_fasta(csv_file, db_fasta)
    txt_to_fasta(txt_file, query_fasta)

    # 创建 BLAST 数据库
    os.system(f"makeblastdb -in {db_fasta} -dbtype prot -out protein_db")

    # 运行 BLAST
    run_blast(query_fasta, "protein_db", output_file)

    # 解析 BLAST 输出并保存相似性矩阵
    similarity_matrix = parse_blast_output(output_file, query_fasta, db_fasta)
    
    # 如果矩阵不是空的，保存为 CSV 文件
    if not similarity_matrix.empty:
        similarity_matrix.to_csv(matrix_output_file)
        print(f"相似性矩阵已保存为 {matrix_output_file}")
    else:
        print("BLAST 输出没有匹配，无法生成相似性矩阵。")

if __name__ == "__main__":
    main()
