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

int main( int argc, char ** argv)
{
	MPI_Init(NULL,NULL);

	MPI_Comm world = MPI_COMM_WORLD;
	int rank, worldSize;
	MPI_Comm_size(world,&worldSize);
	MPI_Comm_rank(world,&rank);
	
	char buf1[255];
	char buf2[255];
	memset(buf1,'\0',255);
	memset(buf2,'\0',255);
	
	FILE *fp1, *fp2;
	if ((fp1 = fopen("./shadow.txt", "r")) < 0)
	{
		printf("Error opening file!: shadow.txt\n");
		MPI_Finalize();
		return 1;
	}
	if ((fp2 = fopen("./words.txt", "r")) < 0)
	{
		printf("Error opening file!: words.txt\n");
		MPI_Finalize();
		return 1;
	}

	
	char *user, *id, *salt, *hash;
	char pre[255], suf[255], temp[255], saltid[255];
	memset(pre,'\0',255);
	memset(suf,'\0',255);
	memset(temp,'\0',255);	
	memset(saltid, '\0', 255);

	int counts = 235888/worldSize;
	int extra = 0;
	if (rank == worldSize - 1)
	{
		extra = 235888 % worldSize;
	}
	printf("Counts = %d\n\n", counts);
	int reads = 1;
	while (reads != rank*counts && rank != 0)
	{
		fscanf(fp2, "%s", buf2);
		reads++;
	}
	char num4[5], num3[4], num2[3];
	memset(num4, '0', 4);
	num4[4] = '\0';
	memset(num3, '0', 3);
	num3[3] = '\0';
	memset(num2, '0', 2);
	num2[2] = '\0';
	
	
	//TODO::::
	//Prettier Printouts and Final Considerations!!!
	// 1. Each word takes forever to try 11111 combinations of numbers
	//    and even longer to perform the sprintf operations to get
	//    each guess.
	// 2. We want to keep an array on the root node 0 that is the length
	//    of the shadow file (in this case, 11). When a password has been
	//    found, it's corresponding position will be updated to a 1, and
	//    it will automatically be 0, for password not found.
	//    - We may also want to keep an array of strings that denote the 
	//      actual password.
	// 3. Need to perform an MPI_Receive, Gather, etc. to bring the results
	//    back to the root. 
	// 4. Need to perform some MPI_File_open to make an output file,
	//    use MPI_DataTypes, MPI_Type_vector, MPI_Type_commit, MPI_Type_create_resized
	//    MPI_File_set_view, and MPI_File_write.
	// 5. Once the formatting is 'okay', time to do some basic timing. Not a happy
	//    moment, considering I am currently running this and it's taking forever...
	//    - Maybe just time the time it takes to complete 1 word, 2 words, depending
	//      on the number of nodes. And estimate the final compile time if we were
	//      to use 10 nodes (since it will take forever). Then we can just take the 
	//      final time and divide by worldSize to see if we are achieving a complete
	//      theoretical speedup.
	// 6. ReadMe questions and such..
	// 
	// (maybe for better readability)
	// int count = 0
	// if (count % 100 == 0)
	//	file_write("Rank x is on 100th word: word")
	//
	// if (guess = password)
	// 	file_write("Password for user is guess")
	// ...
	// if (no password appears..)
	// 	file-write("No password found for user")
	// 	
	// Alright new plan... We're gonna localize the hell out this thing for speed.
	// i.e. remove a ton of fscans, hardcode salts, etc for efficieny	

	for (int j = rank*counts ; j < (rank+1)*counts + extra; j++)
	
		fscanf(fp2, "%s", buf2);
		if (j % 1000 == 0)
		{
			printf("Rank == %d, j == %d, word == %s\n", rank, j, buf2);
		}
		
		for (int k = 0; k < 9999; k++)
		{
			if (k < 1000 && k != 0)
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
			if (k < 100 && k != 0)
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
			if (k < 10 && k != 0)
			{
				if (num2[1] == '9')
				{
					num2[1] == '0';
				}
				else
				{
					num2[1]++;
				}
			}
			fseek(fp1, 0, SEEK_SET); 
			
			for (int i = 0; i < 11; i++) 
			{
				fscanf(fp1, "%s", buf1); 
				user = strtok(buf1,":");
				id = strdup(user);
				id = strtok(NULL,"$");
				salt = strdup(id);
				salt = strtok(NULL,"$");
				hash = strdup(salt);
				hash = strtok(NULL,"$");
				
				snprintf(saltid, sizeof(saltid), "$%s$%s$%s", id, salt, hash); 
				snprintf(temp, sizeof(temp), "$%s$%s", id, salt);
				
				sprintf(pre,"%d%s", k, buf2);
				sprintf(suf,"%s%d", buf2, k);	
				
				if (strcmp(crypt(pre, temp), saltid) == 0) 	
				{
					printf("Password for %s is: %s\n", user, pre);
				}
				else if (strcmp(crypt(suf, temp), saltid) == 0)
				{
					printf("Password for %s is: %s\n", user, suf);
				}
				if (k < 1000)
				{
					sprintf(pre, "%s%s", num4, buf2);
					sprintf(suf, "%s%s", buf2, num4);

					if (strcmp(crypt(pre, temp), saltid) == 0)
					{
						printf("Password for %s is: %s\n", user, pre);
					}
					else if (strcmp(crypt(suf, temp), saltid) == 0)
					{
						printf("Password for %s is: %s\n", user, suf);
					}
				}
				if (k < 100)
				{
					sprintf(pre, "%s%s", num3, buf2);
					sprintf(suf, "%s%s", buf2, num3);

					if (strcmp(crypt(pre, temp), saltid) == 0)
					{
						printf("Password for %s is: %s\n", user, pre);
					}
					else if (strcmp(crypt(suf, temp), saltid) == 0)
					{
						printf("Password for %s is: %s\n", user, suf);
					}
				}
				if (k < 10)
				{
					sprintf(pre, "%s%s", num2, buf2);
					sprintf(suf, "%s%s", buf2, num2);

					if (strcmp(crypt(pre, temp), saltid) == 0)
					{
						printf("Password for %s is: %s\n", user, pre);
					}
					else if (strcmp(crypt(suf, temp), saltid) == 0)
					{
						printf("Password for %s is: %s\n", user, suf);
					}
				}
			}
		}	
	}

	
	fclose(fp1);
	fclose(fp2);
	MPI_Finalize();
	return 0;
}
