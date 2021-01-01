# Unix Basics

## Connecting to `state`

You will first need to connect through `acc.ohsu.edu`

```
ssh laderast@acc.ohsu.edu
```

Once you're in there:

```
ssh state
```

## `pwd` - where am I?

When you sign into your account, you'll be in your *home directory* - this is a directory that houses your files

What is the absolute path to here? You can use `pwd` to find it.

## File Permissions

`chmod`

## Finding Software

`which R`

## Environment Variables

How does Unix know where software is?

It has to do with providing an *environment* variable.

```
echo $PATH
```

## Adding software to your `PATH` using a `.bashrc` file

```
export PATH="$HOME/FastQC/:$PATH"
```

## Editing your `.bashrc` using `nano`




## `which`

## `find`

## `tldr`

# [Command Line Environment](https://missing.csail.mit.edu/2020/command-line/)

## Dotfiles

## Job Control: `ps` and `kill`

## Aliases

## Very Helpful: `tmux`

# Dependencies Manager: `conda` and `miniconda`

