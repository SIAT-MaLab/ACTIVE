import duckdb
import os

def main():
    input_file = 'all_vs_all.tsv'
    output_file = 'genome_similarity.abc'
    db_file = 'diamond_analysis.duckdb'
    memory_limit = '128GB' 
    if os.path.exists(db_file):
        os.remove(db_file)
        
    con = duckdb.connect(db_file)
    
    con.execute(f"SET memory_limit='{memory_limit}';")
    con.execute("SET temp_directory='.';")

    con.execute(f"""
        CREATE TABLE raw_data AS 
        SELECT 
            column0 AS q,
            column1 AS s,
            column2 AS bitscore,
            regexp_replace(column0, '_[^_]+$', '') AS g_query,
            regexp_replace(column1, '_[^_]+$', '') AS g_target
        FROM read_csv('{input_file}', delim='\t', header=False, columns={{'column0':'VARCHAR', 'column1':'VARCHAR', 'column2':'DOUBLE'}});
    """)
    
    con.execute("""
        CREATE TABLE best_hits AS
        SELECT 
            q, g_query, g_target, MAX(bitscore) as max_score
        FROM raw_data
        GROUP BY q, g_query, g_target
    """)

    con.execute("DROP TABLE raw_data")

    con.execute("""
        CREATE TABLE genome_interactions AS
        SELECT 
            g_query, g_target, SUM(max_score) as cumulative_score
        FROM best_hits
        GROUP BY g_query, g_target
    """)

    result_query = """
        COPY (
            WITH self_scores AS (
                SELECT g_query, cumulative_score as self_score
                FROM genome_interactions
                WHERE g_query = g_target
            )
            SELECT 
                t1.g_query, 
                t1.g_target, 
                t1.cumulative_score / t2.self_score AS similarity
            FROM genome_interactions t1
            JOIN self_scores t2 ON t1.g_query = t2.g_query
            WHERE t1.cumulative_score / t2.self_score > 0
        ) TO '{}' (DELIMITER '\t', HEADER FALSE);
    """.format(output_file)
    
    con.execute(result_query)

    
    con.close()
    if os.path.exists(db_file):
        os.remove(db_file)

if __name__ == "__main__":
    main()

