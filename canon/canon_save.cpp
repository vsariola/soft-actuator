// canon_save.cpp : Defines the entry point for the console application.
//

#include <stdlib.h>
#include <stdio.h>
#include <windows.h>
#include "EDSDK.h"
#include "EDSDKErrors.h"
#include "EDSDKTypes.h"

//#pragma comment(lib,"EDSDK.lib");


EdsError downloadImage(EdsDirectoryItemRef directoryItem,char* str){
	EdsError err = EDS_ERR_OK;
	EdsStreamRef stream = NULL;
	// Get directory item information
	EdsDirectoryItemInfo dirItemInfo;
	err = EdsGetDirectoryItemInfo(directoryItem, & dirItemInfo);

	char pname[150];
	//strcpy(pname,"D:\\ContiPre\\temp\\");
	strcpy_s(pname,str);
	strcat_s(pname,dirItemInfo.szFileName);
	//printf(pname);
	
	// Create file stream for transfer destination
	if(err == EDS_ERR_OK) {
		err = EdsCreateFileStream(pname,kEdsFileCreateDisposition_CreateAlways,kEdsAccess_ReadWrite, &stream);
	}
	// Download image
	if(err == EDS_ERR_OK) {
		err = EdsDownload( directoryItem, dirItemInfo.size, stream);
	}
	// Issue notification that download is complete
	if(err == EDS_ERR_OK) {
		err = EdsDownloadComplete(directoryItem);
	}
	// Release stream
	if( stream != NULL) {
		EdsRelease(stream);
		stream = NULL;
	}
	return err;
}

//int _tmain(int argc, _TCHAR* argv[])
int main(int argc,char* argv[])
{
	if ( argc != 3 ){
		printf("two inputs, thank you! %i\n", argc);
		return 0;
	}
	// count input
	EdsUInt32 count;
	count = atoi(argv[1]);
	//count = _ttoi(argv[1]);
	printf("count: %i\n",count);

	char *str = argv[2];
	//CStringA sAnsiString(argv[2]);
	//char *str = sAnsiString;

	//EdsUInt32 count = arg
	//char path;
	
	
	EdsError err = EDS_ERR_OK;
	bool isSDKLoaded = false;
	bool camerafound = false;
	//err = EdsInitializeSDK();
	err = EdsInitializeSDK();
	if(err == EDS_ERR_OK){
		
		isSDKLoaded = true;
		
	} else {
		printf("EDSDK load FAILED \n");
		return 0;
	}
	EdsCameraListRef cameraListRef;
    
	// Get camera list
	err = EdsGetCameraList(&cameraListRef);
	if(err == EDS_ERR_OK){
		camerafound = true;
		
	} else {
		EdsTerminateSDK();
		printf("cameralist load FAILED \n");  
		return 0;
	}


	//the number of connected cameras
    EdsUInt32 numCameras = 0;

	// the number of images in folder
	EdsUInt32 numImages = 0;
    
    //get the number of connected cameras
    err = EdsGetChildCount(cameraListRef, &numCameras);
	//mexPrintf("%i\n",numCameras);
	if(numCameras == 0){
		printf("No camera found \n");
		EdsTerminateSDK();
		return 0;
	}
    
	EdsBaseRef camera;
	err = EdsGetChildAtIndex(cameraListRef, 0, &camera);
	
	err = EdsOpenSession(camera);
	
	EdsVolumeRef volume;
	EdsDirectoryItemRef folder1;
	EdsDirectoryItemRef folder2;
	EdsDirectoryItemRef file;
	EdsGetChildAtIndex(camera, 0 , &volume);
	EdsGetChildAtIndex(volume,0 , &folder1);
	EdsGetChildAtIndex(folder1,0 , &folder2);
	EdsGetChildCount(folder2, &numImages);
	for(unsigned int i = numImages-count; i < numImages; i++){
		err = EdsGetChildAtIndex(folder2, i , &file);
		
		if(err == EDS_ERR_OK){
			downloadImage(file,str);
			EdsDeleteDirectoryItem(file);
		} else {
			break;
		}
	}
	
	
	if(camera != NULL){
		EdsRelease(camera);
	}
	err = EdsCloseSession(camera);
	
    // Terminate SDK
	if(isSDKLoaded)
    {
		EdsTerminateSDK();
	}
	//free(str);
	return 1;
}