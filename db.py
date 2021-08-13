import sqlite3
import pandas


'''
This builds the database from the pickled dataframe
However, it takes a long time (4+ hours)
'''
db = sqlite3.connect("data/db")

#Make the tables if required

create_sample_table = """
    CREATE TABLE IF NOT EXISTS samples(
        uid integer PRIMARY KEY AUTOINCREMENT,
        ena_accession text,
        cogid text,
        sample_date text,
        UNIQUE(ena_accession, cogid, sample_date)
        );
    """
create_aa_mutations_table = """
    CREATE TABLE IF NOT EXISTS aa_mutations(
        uid integer PRIMARY KEY AUTOINCREMENT,
        position integer NOT NULL,
        reference text,
        mutation text,
        UNIQUE(position, reference, mutation)
    );
    """
create_pango_table = """
    CREATE TABLE IF NOT EXISTS pango(
        uid integer PRIMARY KEY AUTOINCREMENT,
        lineage text,
        count integer
    );
    """
create_scorpio_table = """
    CREATE TABLE IF NOT EXISTS scorpio(
        uid integer PRIMARY KEY AUTOINCREMENT,
        lineage text,
        count integer
    );
    """

create_mutations_table = """
    CREATE TABLE IF NOT EXISTS mutations(
        uid integer PRIMARY KEY AUTOINCREMENT,
        sample_uid integer NOT NULL,
        aa_uid integer NOT NULL,
        pango_uid integer,
        scorpio_uid integer,
        FOREIGN KEY (sample_uid) REFERENCES samples (uid),
        FOREIGN KEY (aa_uid) REFERENCES aa_mutations (uid),
        FOREIGN KEY (pango_uid) REFERENCES pango (uid),
        FOREIGN KEY (scorpio_uid) REFERENCES scorpio (uid)
    );
    """

cur = db.cursor()
cur.execute(create_sample_table)
db.commit()
cur = db.cursor()
cur.execute(create_aa_mutations_table)
db.commit()
cur = db.cursor()
cur.execute(create_pango_table)
db.commit()
cur = db.cursor()
cur.execute(create_scorpio_table)
db.commit()
cur = db.cursor()
cur.execute(create_mutations_table)
db.commit()

#Read the dataset
data = pandas.read_pickle("data/MUTATIONS-Spike.pkl.gz")

#Populate the tables
for (index, row) in data.iterrows():
    #Samples table
    cur = db.cursor()
    values = (row["ena_accession"], row["cogid"], row["sample_date"])
    cur.execute("INSERT OR IGNORE INTO samples (ena_accession, cogid, sample_date) VALUES (?, ?, ?);", values)

    #AA Mutations table
    mutation = row["mutation"]
    if "del" in mutation or "ins" in mutation or "indel" in mutation:
        position = int(mutation.split("_")[0])
        mutation = "_".join(mutation.split("_")[1:])
        reference = ""
    else:
        reference = str(mutation[0])
        position = int(mutation[1:-1])
        mutation = mutation[-1]
    cur = db.cursor()
    values = (position, reference, mutation)
    cur.execute("INSERT OR IGNORE INTO aa_mutations (position, reference, mutation) VALUES (?, ?, ?);", values)

    #Lineages

    #Pango
    pango = row["pango_lineage"]
    if not pandas.isnull(pango):
        #There is a pango lineage
        
        #Check if the lineage is already in the table
        cur = db.cursor()
        cur.execute("SELECT * FROM pango WHERE lineage = ?;", (pango, ))
        result = cur.fetchone()
        if result is None:
            #This lineage is not in the table so add it
            cur = db.cursor()
            cur.execute("INSERT INTO pango (lineage, count) VALUES (?, ?);", (pango, 0))
        else:
            #Lineage is in the table so update the count
            count = result[2] + 1
            uid = result[0]
            cur = db.cursor()
            cur.execute("UPDATE pango SET count = ? WHERE uid = ?", (count, uid))

    #Scorpio
    scorpio = row["scorpio_call"]
    if not pandas.isnull(scorpio):
        #There is a scorpio lineage
        
        #Check if the lineage is already in the table
        cur = db.cursor()
        cur.execute("SELECT * FROM scorpio WHERE lineage = ?;", (scorpio, ))
        result = cur.fetchone()
        if result is None:
            #This lineage is not in the table so add it
            cur = db.cursor()
            cur.execute("INSERT INTO scorpio (lineage, count) VALUES (?, ?);", (scorpio, 0))
        else:
            #Lineage is in the table so update the count
            count = result[2] + 1
            uid = result[0]
            cur = db.cursor()
            cur.execute("UPDATE scorpio SET count = ? WHERE uid = ?", (count, uid))
    
    db.commit()
    print("inital: ", index)

for (index, row) in data.iterrows():
    #Get the sample UID
    values = (row["ena_accession"], row["cogid"], row["sample_date"])
    cur = db.cursor()
    cur.execute("SELECT uid FROM samples WHERE ena_accession = ? AND cogid = ? AND sample_date = ?", values)
    result = cur.fetchone()
    if result is None:
        print("????? sample", index, row)
        break
    else:
        sample_uid = result[0]

    #AA mutation
    mutation = row["mutation"]
    if "del" in mutation or "ins" in mutation or "indel" in mutation:
        position = int(mutation.split("_")[0])
        mutation = "_".join(mutation.split("_")[1:])
        reference = ""
    else:
        reference = str(mutation[0])
        position = int(mutation[1:-1])
        mutation = mutation[-1]
    cur = db.cursor()
    values = (position, reference, mutation)
    cur.execute("SELECT uid FROM aa_mutations WHERE position = ? AND reference = ? AND mutation = ?;", values)
    result = cur.fetchone()
    if result is None:
        print("????? aa mutation", index, row)
        break
    else:
        aa_uid = result[0]
    
    #Pango
    pango = row["pango_lineage"]
    if not pandas.isnull(pango):
        cur = db.cursor()
        cur.execute("SELECT uid FROM pango WHERE lineage = ?", (pango, ))
        result = cur.fetchone()
        if result is None:
            print("????? pango", index, row)
            break
        else:
            pango_uid = result[0]
    else:
        pango_uid = None
    #Scorpio
    scorpio = row["scorpio_call"]
    if not pandas.isnull(scorpio):
        cur = db.cursor()
        cur.execute("SELECT uid FROM scorpio WHERE lineage = ?", (scorpio, ))
        result = cur.fetchone()
        if result is None:
            print("????? scorpio", index, row)
            break
        else:
            scorpio_uid = result[0]
    else:
        scorpio_uid = None
    
    cur = db.cursor()
    cur.execute("INSERT INTO mutations (sample_uid, aa_uid, pango_uid, scorpio_uid) VALUES (?, ?, ?, ?);", (sample_uid, aa_uid, pango_uid, scorpio_uid))
    db.commit()
    print("final", index)
    
