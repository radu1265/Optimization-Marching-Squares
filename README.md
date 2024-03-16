# Optimization-Marching-Squares

Algoritmul Marching Squares este paralelizat folosind threads.

Alocarea de memorei a fost mutata din functiile ajutatoare in functia
main. Eliberarea de memorie se face tot in main.

Variabilele sunt trimise in thread_function prin doua structuri de date,
una care pastreaza id-ul thread-ului si un pointer catre o variabila de tip thread_data.
thread_data este celalta structura care mentine toate informatiile necesare pentru rularea functiei,
dar care sunt independente de thread-uri.

thread_function apeleaza functiile paralelizate si utilizeaza pthread_barrier pentru a asigura
sincronizarea thread-urilor.

Functiile paralelizate sunt:rescale_image, sample_grid si march.
Procesul de paralelizare consta in impartirea sarciniilor in cate thread-uri cere utilizatorul.
In fiecare functie paralelizata for-urile care permit aceasta impartire sunt distribuite pe cate un thread.
Calculul de inceput si final al for-ului este acela din labrator.
