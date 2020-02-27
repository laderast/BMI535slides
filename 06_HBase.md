# What is HBase?

- Column store database
- Created as part of Apache project

# What is a column store database?



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

You make a column family by thinking about what columns should go together.

Example:

Customer Table

- Customer ID

Customer Info Family
key: Customer Id
- Customer Name
- Age
- Honorific

Address Column Family:
key: Zipcode
- Address
- City
- State

# Column Families and Relational Databases

- Column families are analogous to single tables in a relational database structure
- Each row key and tuple is analogous to a *row* in a relational database

# The Row Key is everything

The row key determines the position of a row in the database.

Pick it carefully according to the queries you want to make.

# Example of a Row Key

Say we want to 

Table of websites and information
    - row key: domain 

# What's crazy about HBase

- The row key can be (almost) anything: 
    - alphanumeric, integer, even other data structures
    

# 

# Setup

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

# Listing what exists

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

# Creating an HBase Table

```
create 'laderast:gtable`
```

# Describe

```
describe 'laderast:gtable'
```

# HBase versus Hive

# Helpful Links

- [Understanding HBase and BigTable](https://dzone.com/articles/understanding-hbase-and-bigtab)
- [HBase Reference Guide](http://hbase.apache.org/book.html), especially [HBase Data Model](http://hbase.apache.org/book.html#datamodel)
- [Why column stores?](https://blog.pythian.com/why-column-stores/)
- [The beauty of column oriented data](https://towardsdatascience.com/the-beauty-of-column-oriented-data-2945c0c9f560)
