#! /bin/bash

# ECHO COMMAND
echo Hello world!

# VARIABLES: Uppercase by convention, allow letters, numbers, underscores
# Remember to not put space before and after "="
NAME="Wanying"
echo
echo "# Variable demo:"
echo $NAME

# USER INPUT
echo
echo "# User input demo:"
read -p "Enter your name: " NAME
echo "My name is $NAME"

# IF STATEMENT
# [] is a program and need space between the rest (arguments)
# Use = instead of ==
echo
echo "# If statement demo:"
if [ "$NAME" = "Wanying" ]
then
	echo "My name is Wanying"
elif [ "$NAME" = "wanying" ]
then
	echo "My name is wanying"
else
	echo "My name is not Wanying or wanying"
fi

# FOR LOOP
echo
echo "# For-loop demo:"
NAMES="Wanying wanying"
for NAME in $NAMES
do
	echo "My name is $NAME"
done
# Another way of for loop, 
# {start..end..step} is similar as range() function in python
echo "# Another for loop demo, use {start..end..step}"
NAME="Wanying"
for i in {1..10..2}
do
	echo "The number is" $i
	# Use ${var_name} with string defined by "" to do concatenation
	# Sinlge quote will not evaluate variable name in the {}
	echo "The name is ${NAME}, number is ${i}"
done


# WHILE LOOP (read from a file line by line)
# LINE is the line number, CURRENT_LINE is the content of the line
# test.txt is the file to be read in
# Code starts from the next line
# LINE=1
# while read -r CURRENT_LINE
#	do
#		echo "$LINE: $CURRENT_LINE"
#		((LINE++))
# done < "./test.txt"

# FUNCTION
echo
echo "# Function demo:"
greet() {
	echo "First parameter $1, second parameter $2"
}
greet hello world
