# An Introduction to Column Store Databases

## What is HBase?

- Column store database
- Created as part of Apache project

## What is a column store database?

A database technology that focuses on the columns of the data, not the rows.

## Glossary of Terms

From https://www.bmc.com/blogs/hadoop-hbase/

<table cellpadding="1">
<tbody>
<tr>
<td>Row Key</td>
<td>TimeStamp</td>
<td colspan="2">Column Family</td>
<td colspan="2">Column Family</td>
</tr>
<tr>
<td></td>
<td></td>
<td>Column 1</td>
<td>Column 2</td>
<td>Column 3</td>
<td>Column 4</td>
</tr>
</tbody>
</table>

- **Table** - a collection of rows.
- **Row** - a collection of column families.
- **Column Family** - a collection of columns. HBase stores data by column family and not row. This makes for faster retrieval of columns since they are located near each other.
- **Column** - individual data items, like product price. 
- **Timestamp** HBase indexes row keys, columns, and timestamp. The timestamp can be any time format including simple integers which are not time at all. Because timestamp is part of the key, you can have duplicate rows. That means you can have multiple versions of a row. For example you could have an accounting transaction at time n and the same information at some other time.
- **Row Key** - The row key is whatever you store in the row key column. So it’s not just an ordinal number. 
- **Cell** - Google says that HBase is a “sparse, consistent, distributed, multidimensional, sorted map.” Data is stored as this map ((rowkey, column family, column, timestamp) -> value). So you can say that a cell is a column value. That’s not the exact technical definition but an easy way to think about it.

## Why use HBase?

- Optimized for Querying
    - Depending on the data, can be orders of magnitude faster than relational data
- Binary search over data stored in a column is very fast

## Why use HBase?

- Your users are distributed over a large network 
    - HDFS allows for distributed (sharded) storage
    - load balancing of access
- Your data is huge (millions -> billions of rows)
- Your data can be sparse (not all rows will have the same column entries)
- Data is sorted by row key

## Number one use case

Facebook messenger: https://www.facebook.com/notes/facebook-engineering/the-underlying-technology-of-messages/454991608919/

Allows for fast retrieval of messages by timestamp

## Sparsity

Unlike relational databases, NULL values do not take up space.

## How is the Data Stored?

<table cellpadding="1">
<tbody>
<tr>
<td>Row Key</td>
<td>TimeStamp</td>
<td colspan="2">Column Family</td>
<td colspan="2">Column Family</td>
</tr>
<tr>
<td></td>
<td></td>
<td>Column 1</td>
<td>Column 2</td>
<td>Column 3</td>
<td>Column 4</td>
</tr>
</tbody>
</table>

- Data is *sorted and stored* as a sorted data structure (a map)
    - To get the value of a cell, you need this information
        (rowkey, column family, column, timestamp)
- Sorting is critical for fast search and access
- HBase is made for distributed file systems, in particular, Hadoop.

## Your Turn

Identify the row key, the column families, and the columns in the following table:

<table cellpadding="1">
<tbody>
<tr>
<td>Transcript Name</td>
<td>TimeStamp</td>
<td colspan="3">Gene Information</td>
<td colspan="5">Transcript</td>
</tr>
<tr>
<td></td>
<td></td>
<td>Gene ID</td>
<td>GeneSym</td>
<td>Biotype</td>
<td>Transcript ID</td>
<td>Chr</td>
<td>Start</td>
<td>End</td>
<td>Start</td>
<td>Strand</td>
</tr>
</tbody>
</table>

## HBase is a little primitive

Compared to Relational Databases, we often use *denormalized* data because it is faster to access.

## Versioning

A *cell* not only contains a value, but also a *timestamp*. 

Depending on the number of previous versions you store, you may be able to roll back a cell.

## Column Families

An HBase table is decomposed into *column families*, which you can think of as subtables within a larger table.

A *column family* is defined as a key-value pair: *row key* and associated columns

A column family represents the basic unit for adding data. 

## How to design column families

You make a column family by thinking about what columns should go together to speed up your query.

## Example

- row key: chromosome-start position

## Column Families and Relational Databases

- Column families are analogous to single tables in a relational database structure
- Each row key and tuple in a column family is analogous to a *row* in a relational database

## The Row Key is everything

The row key determines the position of a row in the database.

Pick it carefully according to the queries you want to make.

## What's crazy about HBase

- The row key can be (almost) anything: 
    - alphanumeric, integer, even other data structures
- Again, choose it carefully, because it determines performance 

## Setup on `state`

Add these lines to your `~/.bashrc`:

```
#this is where the hbase shell lives
export HBASE_HOME="/home/courses/BMI535/students/hbase/hbase-1.2.6/"  
export PATH="$PATH:$HBASE_HOME/bin"
#need to add JAVA_HOME to make sure it runs
export JAVA_HOME=/usr/lib/jvm/java-8-openjdk-amd64/
export PATH="$PATH:$JAVA_HOME/bin"
```

Remember to

```
source ~/.bashrc
```

## Opening up the HBase Shell

Once you have the above setup, you can open the HBase Shell by running:

```
hbase shell
```

## `list`ing what exists

```
list
```

## Creating your own namespace

A *namespace* is a place of your very own in the HBase database. It helps you avoid what are called *namespace* collisions.

One example of a *namespace collision* is if multiple users tried to create a table called `Gene` in the main workspace. There would be lots of errors, especially if the table structure was different.

So you can create your own namespace in HBase by running

```
create_namespace laderast
```

So, create your own namespace (use your username)!

## Creating an HBase Table / Column Family

If the table doesn't exist, you have to create both at once.

```
create 'laderast:Transcript', 'Gene Information', 'Transcript Information'
```

# What did we just do?

```
describe 'laderast:Transcript'
```

## Data Manipulation Verbs in HBase

- `put` - Insert Data into Table by row key
- `get` - retrieve row by row key
- `scan` - search by row_key or field



## `get`ing data from a table



```
get <'tablename'>, <'rowname'>, {< Additional parameters>}
```

## `scan`ning a table

```
scan <table>, {attributes => ‘value’}
```
## Useful attributes for `scan`

COLUMNS => 'personal_data:name', 
LIMIT => 10, 
STARTROW => '3'
TIMESTAMP => 


## 'put'ing data into `laderast:Transcript`


```
put <'tablename'>,<'rowname'>,<'column_family:column'>,<'value'>
```
```
put 'laderast:Transcript', 'TSPAN6-201' , 'Gene Information: GeneId', 'ENSG00000000003'
put 'laderast:Transcript', 'TSPAN6-201' , 'Gene Information: GeneName', 'TSPAN6'
put 'laderast:Transcript', 'TSPAN6-201' , 'Transcript Information: TranscriptID', 'ENST00000373000'
```

Check that we did something:

```
scan 'laderast:Transcript'
```

## Your turn

Add the following data into your own `laderast:Transcript` table:

In bash (you can use `exit` to get out of `hbase shell`):

1. look in `/home/courses/BMI535/students/hbase/transcript.txt` - copy this file to your own directory
2. You'll need to modify the `Transcript` table name to fit your own namespace. Hint: search and replace.
3. When you're done, try running: 

```
hbase shell transcript.txt
```


# Table Design Rules of Thumb

[From HBase Reference](http://hbase.apache.org/book.html): 

- Aim to have regions sized between 10 and 50 GB.

- Aim to have cells no larger than 10 MB, or 50 MB if you use mob. Otherwise, consider storing your cell data in HDFS and store a pointer to the data in HBase.

- A typical schema has between 1 and 3 column families per table. HBase tables should not be designed to mimic RDBMS tables.

- Around 50-100 regions is a good number for a table with 1 or 2 column families. Remember that a region is a contiguous segment of a column family.

- Keep your column family names as short as possible. The column family names are stored for every value (ignoring prefix encoding). They should not be self-documenting and descriptive like in a typical RDBMS.

- If you are storing time-based machine data or logging information, and the row key is based on device ID or service ID plus time, you can end up with a pattern where older data regions never have additional writes beyond a certain age. In this type of situation, you end up with a small number of active regions and a large number of older regions which have no new writes. For these situations, you can tolerate a larger number of regions because your resource consumption is driven by the active regions only.

- If only one column family is busy with writes, only that column family accomulates memory. Be aware of write patterns when allocating resources.

# Column Family Tips

[From HBase Reference](http://hbase.apache.org/book.html):

- HBase currently does not do well with anything above two or three column families so keep the number of column families in your schema low.

# HBase versus Hive

- https://www.xplenty.com/blog/hive-vs-hbase/

# HBase in Bioinformatics and Clinical Informatics

- https://www.nitrc.org/forum/message.php?msg_id=21408
= https://www.hindawi.com/journals/cmmm/2017/6120820/


# Helpful Links

- [Understanding HBase and BigTable](https://dzone.com/articles/understanding-hbase-and-bigtab)
- [Apache HBase: Overview and Use Cases](https://events.static.linuxfound.org/sites/events/files/slides/ApacheBigData2016.pdf)


- https://www.slideshare.net/HBaseCon/case-studies-session-7
- [HBase Beginners Guide](https://acadgild.com/blog/hbase-tutorial-beginners-guide)
- [HBase Reference Guide](http://hbase.apache.org/book.html), especially [HBase Data Model](http://hbase.apache.org/book.html#datamodel)
- [Why column stores?](https://blog.pythian.com/why-column-stores/)
- [The beauty of column oriented data](https://towardsdatascience.com/the-beauty-of-column-oriented-data-2945c0c9f560)
