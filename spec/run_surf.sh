#!/bin/bash

# Putanja do vašeg izvršnog programa
EXECUTABLE="./surf.ln"

# Pitajte korisnika za naziv slike (bez ".jpg" ekstenzije)
echo "Unesite naziv slike iz data foldera(bez .jpg ekstenzije):"
read IMAGE_NAME

# Argumenti koje želite proslediti programu
ARGS="-i ${IMAGE_NAME}.jpg -o ../data/${IMAGE_NAME}.surf"

# Pokretanje programa
$EXECUTABLE $ARGS

# Dodajte sve ostale argumente koje želite proslediti programu
# ARGS="$ARGS -thress 1000 -d -ms 3 -oc 4 -ss 2 -u -e -in 4 -p1 regions.txt"
# $EXECUTABLE $ARGS

