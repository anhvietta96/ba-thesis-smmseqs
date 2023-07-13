Die Berechnung in multisets_encode basiert im der Summierung verschiedener
Gewichte aus einer Tabelle, die durch ein Python-Script berechnet wird. 
Ich weiss die Details nicht mehr. Daher ist Ihre Aufgabe hier auch, die 
Funktionsweise herauszufinden und sie spaeter in der B.Sc. zu beschreiben.
Um den Pythoncode zu verstehen, m"ussen Sie vermutlich auch die PfN1-Folien 
zu den Multisets lesen.

Vorgehen:

- make 

  aufrufen, um multisets.o zu kompilieren.

- Ein C-hauptprogramm schreiben, dass alle Multisets fuer ein Alphabet und ein 
  q aufzaehlt, jedes mit der encode funktion aus multisets.c codiert, die 
  Codes aufsammeln und sortieren. Dann muss sich die Folge der Code bis zum 
  maximalen Code ergeben.

- Eine C++-Klasse implementieren, die als template Parameter die q-gram 
  laenge <= 7 hat und mit Alphabetgroesse = 20 funktioniert, so dass sie
  die gleichen Tabellen verwenden k"pnnen.

- Die C++-Klasse um eine Template-Variable f"ur die ALphabetgr"osse erweitern
  und die Berechnungen der Gewichte aus den Python-Scripten durch C++-Code
  in der Klasse ersetzen.

VG. Stefan Kurtz
