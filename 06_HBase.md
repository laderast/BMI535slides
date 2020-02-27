# Why use HBase?

- Storage and writing is very fast

# Column Families

# Setup

Add to your `~/.bashrc`:

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

# Creating your own namespace



```
