//Author: Mitchell Krystofiak
//Class: COSC 420 - Project 1: Password Cracker
//Date: October 17, 2021
//Description: Given unique hashes, attemp to decipher the salt and id of the 
//	       alogirithm used to create an encrypted password.

#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>
#include<time.h>
#include<string.h>
#include<unistd.h>
#include<crypt.h>
#include<stdbool.h>

void massfree(char **username, char *ids, char **salts, char **hashes, char *buf2)
{
	free(buf2);
	for (int i = 0; i < 11; i++)
	{
		free(username[i]);
		//free(ids[i]);
		free(salts[i]);
		free(hashes[i]);
	}
	free(username);
	free(ids);
	free(salts);
	free(hashes);
}

int main( int argc, char ** argv)
{
	//MPI INITIALIZATION AND VARIABLES
	//--------------------------------

	MPI_Init(NULL,NULL);

	MPI_Comm world = MPI_COMM_WORLD;
	int rank, worldSize;
	MPI_Comm_size(world,&worldSize);
	MPI_Comm_rank(world,&rank);
	double t1, t2, t3, t4;	
	
	//INITIALIZE ARRAYS FOR PASSWORD CRACKING
	//---------------------------------------
	
	t1 = MPI_Wtime();
	char **username = malloc(sizeof(char*)*11);
	char *ids      = malloc(sizeof(char)*11);
	char **salts    = malloc(sizeof(char*)*11);
	char **hashes   = malloc(sizeof(char*)*11);
	for (int i = 0; i < 11; i++)
	{
		username[i] = malloc(sizeof(char)*20); //max length plus extra
		//ids[i]      = malloc(sizeof(char)*1);  //1
		salts[i]    = malloc(sizeof(char)*5);  //'ab'
		hashes[i]   = malloc(sizeof(char)*40); //say, pre+suf+32 length words = 40
	}

	//LOAD 'shadow.txt' AND RIP EACH LINE APART
	//-----------------------------------------

	FILE *fp1;
	if ((fp1 = fopen("./shadow.txt", "r")) == NULL)
	{
		printf("Error opening file!: shadow.txt\n");
		MPI_Finalize();
		return 1;
	}
	
	char* buf1 = malloc(100*sizeof(char)); 
	char *user, *id, *salt, *hash;
	int  user_num = 0;

	while ((fscanf(fp1, "%s", buf1)) > 0)
	{
		user = strtok(buf1,":");
		strcpy(username[user_num], user); 
		id = strdup(user);
		id = strtok(NULL,"$");
		strcpy(&ids[user_num], id);
		salt = strdup(id);
		salt = strtok(NULL, "$");
		strcpy(salts[user_num], salt);
		hash = strdup(salt);
		hash = strtok(NULL, "\n");
		snprintf(hashes[user_num], 40, "$%s$%s$%s", id, salt, hash);
		user_num++;
	}
	fclose(fp1);
	free(buf1);

	//LOAD 'words.txt' AND PREPARE COUNTS
	//-----------------------------------
	
	FILE *fp2;
	if ((fp2 = fopen("./words.txt", "r")) == NULL)
	{
		printf("Error opening file!: words.txt\n");
		MPI_Finalize();
		return 1;
	}
	
	//LOCALIZE ALL DICTIONARYS FOR SPEED UP
	char* buf2 = malloc(sizeof(char)*100);	
	char pre[50], suf[50], pre1[50], suf1[50], pre2[50], suf2[50], pre3[50], suf3[50];
	memset(pre, '\0', 50);
	memset(suf, '\0', 50);
	memset(pre1, '\0', 50);
	memset(suf1, '\0', 50);
	memset(pre2, '\0', 50);
	memset(suf2, '\0', 50);
	memset(pre3, '\0', 50);
	memset(suf3, '\0', 50);

	char num4[5], num3[4], num2[3];
	memset(num4, '0', 4);
	num4[4] = '\0';
	memset(num3, '0', 3);
	num3[3] = '\0';
	memset(num2, '0', 2);
	num2[2] = '\0';

	int counts = 235888/worldSize;
	int extra = 0;

	if (rank == worldSize - 1)
	{
		extra = 235888 % worldSize;
	}
	char * Ldict = malloc(sizeof(char)*(counts+extra));

	int reads = 1;
	while (reads != rank*counts && rank != 0)
	{	
		fscanf(fp2, "%s", buf2);
		reads++;
	}
	char saltid[] = "$1$ab";	
	
	//PREPARE FILES FOR READABLE OUTPUT AND LOGS
	//------------------------------------------
	
	FILE  *output;
	//make new file for each
	//change position of sprintfs to before 11, make new arrays 
	if ((output = fopen("./output.txt", "w+")) == NULL)
	{
		printf("Error opening file!: output.txt\n");
		MPI_Finalize();
		return 1;
	}

	//PREPARE ASYNCHRONOUS MESSAGES FOR COMPLETED PASSWORDS
	//-----------------------------------------------------
	
	int *final = malloc(sizeof(int)*11);
	memset(final, 0, sizeof(int)*11);
	
	int allDone = 0;
	int one = 1;
	int done;
	MPI_Request sendRqsts[worldSize*11];
	MPI_Request recvRqsts[worldSize*11];

	//prepares each node to recv incoming messages
	for (int i = 0; i < worldSize; i++)
	{
		for (int j = 0; j < 11; j++)
		{
			MPI_Irecv(&final[j], 1, MPI_INT, i, j,
				world, &recvRqsts[i*11 + j]);
		}
	}

	//BEGIN PASSWORD CRACKING -- SLOW :(
	//----------------------------------

	bool forceStop = false;	

	int totalCounts = 0;
	int count = 0;
	t4 = MPI_Wtime();
	for (int j = rank*counts ; j < (rank+1)*counts + extra && !allDone && !forceStop; j++)
	{
		allDone = 1;
		for (int i = 0; i < 11; i++)
		{
			if ( !final[i]) 
			{
				allDone = 0;
			}
		}
		
		t1 = MPI_Wtime();
		fscanf(fp2, "%s", buf2);
		if (totalCounts % 10 == 0)
		{
			//sprintf(log, "Rank %d is working on word %s\n", rank, buf2);
		}
		totalCounts++;
		for (int k = 0; k < 9999 && !allDone && !forceStop; k++)
		{
			printf("Rank %d total rate: %0.2f\n", rank, ((float)count)/(MPI_Wtime()-t4));
			if (k < 1001 && k != 0)
			{
				if (num4[3] == '9')
				{
					num4[3] = '0';
					if (num4[2] == '9')
					{
						num4[2] = '0';
						if(num4[1] == '9')
						{
							num4[1] = '0';
						}
						else
						{
							num4[1]++;
						}
					}
					else
					{
						num4[2]++;
					}
				} 
				else 
				{
					num4[3]++;
				}
			}
			if (k < 101 && k != 0)
			{
				if (num3[2] == '9')
				{
					num3[2] = '0';
					if (num3[1] == '9')
					{
						num3[1] = '0';
					}
					else 
					{
						num3[1]++;
					}
				}
				else
				{
					num3[2]++;
				}
			}
			if (k < 11 && k != 0)
			{
				if (num2[1] == '9')
				{
					num2[1] = '0';
				}
				else
				{
					num2[1]++;
				}
			}
			
			sprintf(pre,"%d%s", k, buf2);
			sprintf(suf,"%s%d", buf2, k);	
			sprintf(pre1,"%s%s", num4, buf2);
			sprintf(suf1,"%s%s", buf2, num4);	
			sprintf(pre2,"%s%s", num3, buf2);
			sprintf(suf2,"%s%s", buf2, num3);
			sprintf(pre3,"%s%s", num2, buf2);
			sprintf(suf3,"%s%s", buf2, num2);

			for (int i = 0; i < 11 && !allDone; i++) 
			{
				if (final[i])
					continue;
					
				count+=2;
				if (strcmp(crypt(pre, saltid), hashes[i]) == 0) 	
				{
					t3 = MPI_Wtime();
					fprintf(output,"Rank: %d\n", rank);	
					fprintf(output,"Password found for %s is: %s\n", username[i], pre);
					fprintf(output,"Time taken: %f seconds\n\n", t3-t4);
					fflush(output);
					final[i] = 1;
					break;
				}
				else if (strcmp(crypt(suf, saltid), hashes[i]) == 0)
				{
					t3 = MPI_Wtime();
					fprintf(output,"Rank: %d\n", rank);	
					fprintf(output,"Password found for %s is: %s\n", username[i], suf);
					fprintf(output,"Time taken: %f seconds\n\n", t3-t4);	
					fflush(output);
					final[i] = 1;
					break;	
				}
				else if (strcmp(crypt(buf2, saltid), hashes[i]) == 0)
				{
					t3 = MPI_Wtime();
					fprintf(output,"Rank: %d\n", rank);	
					fprintf(output,"Password found for %s is: %s\n", username[i], suf);
					fprintf(output,"Time taken: %f seconds\n\n", t3-t4);
					fflush(output);
					final[i] = 1;
					break;	
				}
				if (k < 1000)
				{
					count+=2;
					if (strcmp(crypt(pre1, saltid), hashes[i]) == 0)
					{
						t3 = MPI_Wtime();
						fprintf(output,"Rank: %d\n", rank);
						fprintf(output,"Password found for %s is: %s\n", username[i], pre);
						fprintf(output,"Time taken: %f seconds\n\n", t3-t4);
						fflush(output);
						final[i] = 1;
						break;		
					}
					else if (strcmp(crypt(suf1, saltid), hashes[i]) == 0)
					{
						t3 = MPI_Wtime();
						fprintf(output,"Rank: %d\n", rank);
						fprintf(output,"Password found for %s is: %s\n", username[i], suf);
						fprintf(output,"Time taken: %f seconds\n\n", t3-t4);
						fflush(output);
						final[i] = 1;
						break;	
					}
				}
				if (k < 100)
				{
					count+=2;
					if (strcmp(crypt(pre2, saltid), hashes[i]) == 0)
					{
						t3 = MPI_Wtime();
						fprintf(output,"Rank: %d\n", rank);	
						fprintf(output,"Password found for %s is: %s\n", username[i], pre);
						fprintf(output,"Time taken: %f seconds\n\n", t3-t4);
						fflush(output);
						final[i] = 1;
						break;	
					}
					else if (strcmp(crypt(suf2, saltid), hashes[i]) == 0)
					{
						t3 = MPI_Wtime();
						fprintf(output,"Rank: %d\n", rank);
						fprintf(output,"Password foundfor %s is: %s\n", username[i], suf);
						fprintf(output,"Time taken: %f seconds\n\n", t3-t4);
						fflush(output);
						final[i] = 1;
						break;	
					}
				}
				if (k < 10)
				{
					count+=2;
					if (strcmp(crypt(pre3, saltid), hashes[i]) == 0)
					{
						t3 = MPI_Wtime();
						fprintf(output,"Rank: %d\n", rank);	
						fprintf(output,"Password for %s is: %s\n", username[i], pre);
						fprintf(output,"Time taken: %f seconds\n\n", t3-t4);
						fflush(output);
						final[i] = 1;
						break;		
					}
					else if (strcmp(crypt(suf3, saltid), hashes[i]) == 0)
					{
						t3 = MPI_Wtime();
						fprintf(output,"Rank: %d\n", rank);
						fprintf(output,"Password for %s is: %s\n", username[i], suf);
						fprintf(output,"Time taken: %f seconds\n\n", t3-t4);
						fflush(output);
						final[i] = 1;
						break;
					}
				}
				//Send messages
				for (int m = 0; m < 11; m++)
				{
					if (final[m] == 1)
					{
						for (int n = 0; n < worldSize; n++)
						{
							MPI_Isend(&one, 1, MPI_INT,
								n, m, world, &sendRqsts[n*11 + m]);
						}
					}
				}
				//Test for sent messages
				for (int m = 0; m < 11*worldSize; m++)
				{
					MPI_Test(&recvRqsts[m], &done, MPI_STATUS_IGNORE);
				}
			}//i
		}//k
		t2 = MPI_Wtime();
		//printf("Number of hashes for word %s: %d, time: %f\n", buf2, count, t2 - t1);	
	}//j
	
	MPI_Barrier(world);	
	if (rank == 0)
	{
		fprintf(output, "\nPasswords not discovered:\n");
		for (int i = 0; i < 11; i++)
		{
			if (final[i] == 0)
				fprintf(output, "User: %s, hash: %s\n", username[i], hashes[i]);
		}
		fflush(output);
	}

	fclose(fp2);
	fclose(output);
	massfree(username, ids, salts, hashes, buf2);
	free(final);
	MPI_Finalize();
	return 0;
}
