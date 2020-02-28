# What is HBase?

- Column store database
- Created as part of Apache project

# What is a column store database?



# Glossary of Terms

From https://www.bmc.com/blogs/hadoop-hbase/

- **Table** a collection of rows.
- **Row** a collection of column families.
- **Column Family** a collection of columns. HBase stores data by column family and not row. This makes for faster retrieval of columns since they are located near each other.
- **Column** individual data items, like product price. For example, in a system designed to store product information you could have a column family called characteristics and then another called inventory. Characteristics could contain the columns description, manufacturer, and serial number. Inventory could include itemCount, SKU (stock keeping unit), EAN (barcode Europe), and UPC (barcode USA). You reference the columns like this characteristics:itemCount.
- **Timestamp** HBase indexes row keys, columns, and timestamp. The timestamp can be any time format including simple integers which are not time at all. Because timestamp is part of the key, you can have duplicate rows. That means you can have multiple versions of a row. For example you could have an accounting transaction at time n and the same information at some other time.
- **Row Key** The row key is whatever you store in the row key column. So it’s not just an ordinal number. In this example, we use product number. You do not give the row key a name, like productNumber, like you do with columns. HBase keeps row keys in alphabetical order so that similar items are kept next to each other. For example, in the Google documentation that explains Google Big Table they use the example of their web index. Data for abc.com is kept next to sales.abc.com. (In other to store those next to each other Google stores those domain names backwards like com.abc and com.abc.sales.)
- **Cell** Technically HBase stores data in maps. Python, Scala, and Java programmers know that a map is a (key->value) data structure, with no duplicate keys allowed. Google says that HBase is a “sparse, consistent, distributed, multidimensional, sorted map.” Data is stored as this map ((rowkey, column family, column, timestamp) -> value). So you can say that a cell is a column value.That’s not the exact technical definition but an easy way to think about it.

# Why use HBase?

- Optimized for Querying
    - Depending on the data, can be orders of magnitude faster than relational data
- Binary search over data stored in a column is very fast

# Why use HBase?

- Your data is huge (millions -> billions of rows)
- Data is versioned
- Your data can be sparse (not all rows will have the same column entries)
- Missing data doesn't take memory
- Data is sorted by row key

# Number one use case

Facebook timeline 

# Sparsity

Unlike relational databases, NULL values do not take up space.

# How is the Data Loaded?

HBase is not optimized for constant table updates.

Data is usually bulk loaded.

# How is the Data Stored?

- Data is *sorted and stored* as a sorted data structure
- Sorting is critical for fast search and access
- HBase is made for distributed file systems, in particular, Hadoop.

# HBase is a little primitive

Compared to Relational Databases, we often use *denormalized* data because it is faster to access.

# Versioning

A *cell* not only contains a value, but also a *timestamp*. 

Tables in HBase are versioned.

# The architecture of HBase



# Column Families

A table is decomposed into *column families*, which you can think of as subtables within a larger table.

A column family is defined as a key-value pair: *row key* and associated columns

A column family represents the basic unit for adding data. 

You make a column family by thinking about what columns should go together to speed up your query.

# Example

- row key: chromosome-start position

We have a customer database where we need to quickly identify everyone in a zip code.

## Customer Table
- row key: zipcode-address

### Customer Info Column Family:
- Customer ID
- Customer Name
- Age
- Honorific

### Address Column Family:
- Address
- City
- State

# Column Families and Relational Databases

- Column families are analogous to single tables in a relational database structure
- Each row key and tuple in a column family is analogous to a *row* in a relational database

# The Row Key is everything

The row key determines the position of a row in the database.

Pick it carefully according to the queries you want to make.

# Example of a Row Key

Why is `zipcode-address` used in our example?

# What's crazy about HBase

- The row key can be (almost) anything: 
    - alphanumeric, integer, even other data structures
- Again, choose it carefully, because it determines performance

# Setup on `state`

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

# Opening up the HBase Shell

Once you have the above setup, you can open the HBase Shell by running:

```
hbase shell
```

# `list`ing what exists

```
list
```

# Creating your own namespace

A namespace is a place of your very own in the HBase database. It helps you avoid what are called *namespace* collisions.

One example of a *namespace collision* is if multiple users tried to create a table called `Gene` in the main workspace. There would be lots of errors, especially if the table structure was different.

So you can create your own namespace in HBase by running

```
create_namespace laderast
```

So, create your own namespace (use your username)!

# Creating an HBase Table / Column Family

If the table doesn't exist, you have to create both at once.

```
create 'laderast:gtable', "gene_info"
```

# 'put'ing data into `laderast:gtable`



# Describe

```
describe 'laderast:gtable'
```

# Data Manipulation Verbs in HBase

- put
- get
- scan
- disable
- drop


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

# Helpful Links
- https://www.xplenty.com/blog/hive-vs-hbase/
- https://events.static.linuxfound.org/sites/events/files/slides/ApacheBigData2016.pdf
- [Understanding HBase and BigTable](https://dzone.com/articles/understanding-hbase-and-bigtab)
- [HBase Beginners Guide](https://acadgild.com/blog/hbase-tutorial-beginners-guide)
- [HBase Reference Guide](http://hbase.apache.org/book.html), especially [HBase Data Model](http://hbase.apache.org/book.html#datamodel)
- [Why column stores?](https://blog.pythian.com/why-column-stores/)
- [The beauty of column oriented data](https://towardsdatascience.com/the-beauty-of-column-oriented-data-2945c0c9f560)
