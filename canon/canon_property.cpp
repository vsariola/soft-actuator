#include <stdlib.h>
#include <stdio.h>
#include <windows.h>
#include "EDSDK.h"
#include "EDSDKErrors.h"
#include "EDSDKTypes.h"


//int _tmain(int argc, _TCHAR* argv[])
int main(int argc, CHAR* argv[])
{
	if ( argc != 3 ){
		printf("two input, thank you! %i\n", argc);
		return 1;
	}
	//int c = _ttoi(argv[1]);
	int c = atoi(argv[1]);
	CHAR* TvHex = argv[1];
	CHAR* AvHex = argv[2];

	//int ii = _ttoi(argv[2]);
	//printf("c: %i\n",c);
	EdsError err = EDS_ERR_OK;
	bool isSDKLoaded = false;
	//err = EdsInitializeSDK();
	err = EdsInitializeSDK();
	if(err == EDS_ERR_OK){
		isSDKLoaded = true;
	} else {
		printf("Error\n");
		return 1;
	}
		EdsCameraListRef cameraListRef;
    
	// Get camera list
	err = EdsGetCameraList(&cameraListRef);

	//the number of connected cameras
    EdsUInt32 numCameras = 0;
    
    //get the number of connected cameras
    err = EdsGetChildCount(cameraListRef, &numCameras);
	//mexPrintf("%i\n",numCameras);
    
	EdsBaseRef camera;
	err = EdsGetChildAtIndex(cameraListRef, 0, &camera);
	err = EdsOpenSession(camera);

	EdsDataType dataType;
	EdsUInt32 dataSizeTv = 0;
	EdsUInt32 dataSizeAv = 0;
	EdsVoid * Tv;
	EdsVoid * Av;

	err = EdsSendStatusCommand (camera,kEdsCameraStatusCommand_UILock,0);

	err = EdsGetPropertySize(camera, kEdsPropID_Tv, 0 , &dataType, &dataSizeTv);
	err = EdsGetPropertySize(camera, kEdsPropID_Av, 0 , &dataType, &dataSizeAv);
	if(err == EDS_ERR_OK){
		err = EdsGetPropertyData(camera, kEdsPropID_Tv, 0 , dataSizeTv, &Tv);
		err = EdsGetPropertyData(camera, kEdsPropID_Av, 0 , dataSizeAv, &Av);
	}
	printf("shutterspeed: %x\n",Tv);
	printf("shutterspeed: %x\n",Av);

	if(err == EDS_ERR_OK){
		unsigned int Tvint; 
		unsigned int Avint;
		sscanf(TvHex,"%x", &Tvint);
		sscanf(AvHex,"%x", &Avint);
		EdsUInt32 TvNew  = Tvint; 
		EdsUInt32 AvNew  = Avint; 
		//printf("shutterspeed: %x\n",TvNew);
		err = EdsSetPropertyData(camera, kEdsPropID_Tv, 0 , dataSizeTv, &TvNew);
		//err = EdsGetPropertyData(camera, kEdsPropID_Tv, 0 , dataSizeTv, &Tv);
		err = EdsSetPropertyData(camera, kEdsPropID_Av, 0 , dataSizeAv, &AvNew);
		//err = EdsSetPropertyData(camera, kEdsPropID_Av, 0 , dataSizeAv, &Av);
	}

	//printf("shutterspeed: %x\n",Tv);

	err = EdsSendStatusCommand (camera,kEdsCameraStatusCommand_UIUnLock,0);

	// close session and Release camera
	if(err == EDS_ERR_OK){
		err = EdsCloseSession(camera);
	}
	
	if(camera != NULL){
		EdsRelease(camera);
	}
    // Terminate SDK
	if(isSDKLoaded)
    {
		EdsTerminateSDK();
	}

	return 0;
}
