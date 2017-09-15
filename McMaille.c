/*
*     Version 4.00 parallelized with OpenMP for multi-core processors
*       but this version is slightly modified for monoprocessors
*
*     MAILLE in french = CELL in english 
*     Mc for Monte Carlo
*     Pronounce : MacMy
*
************************************************************************ 
*
*     A Monte Carlo and grid search code for indexing powder patterns
*
*     For more details see the documentation at 
*                   http://www.cristal.org/McMaille/ 
*             or    http://sdpd.univ-lemans.fr/McMaille/ 
*
*              by A. Le Bail - September 2002 for version 0.9
*                              as well as for versions 1.0, 2.0 and 3.0 
*                              October 2006 for version 4.00            
*                        alb@cristal.org 
*                        http://www.cristal.org/ 
*
*                        Résidence Cristal - Appt 213 
*                        2, rue de Gasperi 
*                        72100 Le Mans 
*                        FRANCE
*
*   Versions 0.9 : cubic only
*            1.0 : hexagonal/trigonal/rhombohedral, tetragonal, 
*                  orthorhombic added, plus .ckm and .prf files
*            2.0 : monoclinic and triclinic added in MC
*                  but not in grid search (too long)
*            3.0 : columnar peak shapes instead of Gaussian
*                  in versions 0.9-2.0
*                  no Le Bail fit contrarily to versions 0.9-2.0
*                  only fit by percentage of inclusion of the
*                  calculated column into the observed one
*            3.02: black box mode  
*            3.03: improved Monte Carlo
*            3.04: two-phases mode
*            4.00: automatisation improved : more chances to identify
*                   the correct cell in "black box" mode
*                  Identification of the Bravais lattice
*                  Parallelization by using OpenMP directives
*                  improving the speed with multicore processors
*                  speed x1.7 to 1.8 with "dual core" or "core duo"
*                  speed x3.6 expected with the quad core in 2007
*                  speed x79 expected with the 80-core in 2012...;-)
*                  
*                      
*
************************************************************************ 
*
*    Copyright (C) 2002-2006 Armel Le Bail 
*
* This program is free software; you can redistribute it and/or 
* modify it under the terms of the GNU General Public License 
* as published by the Free Software Foundation.
*
* This program is distributed in the hope that it will be useful, 
* but WITHOUT ANY WARRANTY; without even the implied warranty of 
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License 
* along with this program; if not, write to the Free Software 
* Foundation, Inc.,
* 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*
*
*********************************************************************** 
*
*
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>

#define VERSION "4.00"
#define PROCS 1 //TODO change to actual number of processors
#define PI 114.59156  // 360./3.1415926
#define FILENAME_SIZE 20
#define FILE_WITH_EXTENSION_SIZE FILENAME_SIZE + 4

// Note: This function returns a pointer to a substring of the original string.
// If the given string was allocated dynamically, the caller must not overwrite
// that pointer with the returned value, since the original pointer must be
// deallocated using the same allocator with which it was allocated.  The return
// value must NOT be deallocated using free() etc.
char *trimwhitespace(char *str)
{
	char *end;

	// Trim leading space
	while(isspace((unsigned char)*str)) str++;

	if(*str == 0)  // All spaces?
		return str;

	// Trim trailing space
	end = str + strlen(str) - 1;
	while(end > str && isspace((unsigned char)*end)) end--;

	// Write new null terminator
	*(end+1) = 0;

	return str;
}

FILE *openFile(const char *file_name, const char *ext, const char *mode)
{
	//TODO - change [FILENAME_SIZE + EXTENSION_SIZE]
	char result_file_name[FILE_WITH_EXTENSION_SIZE];
	snprintf(result_file_name, FILE_WITH_EXTENSION_SIZE, "%s%s", file_name, ext);
	FILE *result = fopen(result_file_name, mode);	
	return result;
}

char *getFormattedDate()
{
	time_t timer;
	char *result;
	size_t size = 40;

	result = (char *)malloc(size * sizeof(char));

	struct tm* tm_info;

	time(&timer);
	tm_info = localtime(&timer);

	strftime(result, size, "\n\n %d-%b-%Y\t\t%H hour %M min %S Sec", tm_info);

	return result; 
}

void writeFormattedDate(FILE *file)
{
	char *date = getFormattedDate();
	fwrite(date, strlen(date), 1, file);
	free(date);
}

void writeTotalTime(FILE *file, const time_t *begin, time_t *end)
{
	writeFormattedDate(file);
	time(end);
	double total_time = difftime(*end, *begin);
	char buffer[50];
	snprintf(buffer, sizeof(buffer), "\n\n Total CPU time elapsed in seconds : %.f", total_time);
	fwrite(buffer, strlen(buffer), 1, file);
	printf("\n\nType any character to continue : \n");
	getchar();
}

void insertHeader(FILE *file, const char *file_name)
{
	char buffer[300];
	snprintf(buffer, sizeof(buffer), "\n =================================================================\n McMaille version %s     by A. Le Bail - 2006 - alb@cristal.org\n =================================================================\n\n Using generic filename : %s\n\n   Number of Processors : \t\t%d\n", VERSION, file_name, PROCS);
	fwrite(buffer, strlen(buffer), 1, file);
}


/*
* Open files using name from command line and standard extensions.
* The OPEN statements may need to be changed for some computers.
* Subroutine SXNM gets the generic filename from the command line.
* If nothing is found, then the user is prompted for the filename.
*/
//==============================================================================
char *programInit(int argc, char *argv[])
{
	char *file_name;
	file_name = (char *)malloc(FILENAME_SIZE * sizeof(char));
	
	int write_file_name = (argc != 2);
	do 
	{
		if (write_file_name)
		{
			printf("\nEntry file (no extension) ??");
			scanf("%s", file_name);
			trimwhitespace(file_name);
		}
		else
		{
			snprintf(file_name, FILENAME_SIZE, "%s", argv[1]);
			write_file_name = 1;
		}

		char file_name_ext[FILE_WITH_EXTENSION_SIZE];
		snprintf(file_name_ext, sizeof(file_name_ext), "%s.dat", file_name);
		if ( access(file_name_ext, F_OK) != -1)
		{
			return file_name;
		}

		printf("\n\nThat file does not exist, try again...");
	} while (1);	

}

FILE *createImpFile(char *file_name)
{
	FILE *dat_file = openFile(file_name, ".dat", "r");
	FILE *inp_file = openFile(file_name, ".inp", "w");

	char *buffer = NULL;
	size_t buffsize = 0;
	ssize_t nread;
	while ((nread = getline(&buffer, &buffsize, dat_file)) != -1)
	{
		if (buffer[0] != '!')
		{
			fwrite(buffer, nread, 1, inp_file);
		}
	}

	free(buffer);
	fclose(dat_file);

	printf("\n\nMcMaille version %s\nDataFile : %s\n", VERSION, file_name);
	
	FILE *imp_file = openFile(file_name, ".imp", "w");
	insertHeader(imp_file, file_name);

	fclose(inp_file);

	return imp_file;
}

//==============================================================================

int main(int argc, char *argv[])
{

	//TODO #OMP THREADPRIVATE(/cal/,/cal2/) 

	//Search for the number of processors available 
	//OMP_GET_NUM_PROCS()
	printf("\nNumber of used processors : \t\t%d\n", PROCS);

	int n1 = 1;
	int n2 = 2;
	double x = 0.0;

	char *file_name = programInit(argc, argv);
	FILE *imp_file = createImpFile(file_name);

	time_t time_begin, time_end;

	writeFormattedDate(imp_file);
	time(&time_begin);

	free(file_name);
	fclose(imp_file);

	return 0;
}