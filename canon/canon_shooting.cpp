//#include "stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <windows.h>
#include "EDSDK.h"
#include "EDSDKErrors.h"
#include "EDSDKTypes.h"


//int _tmain(int argc, _TCHAR* argv[])
int main(int argc, CHAR* argv[])
{
	if ( argc != 2 ){
		//printf("one input, thank you! %i\n", argc);
		return 1;
	}
	//int c = _ttoi(argv[1]);
	int c = atoi(argv[1]);

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

	
		if (c == 1){
			//err = EdsSendCommand(camera, kEdsCameraCommand_TakePicture , 0);
			err = EdsSendCommand(camera, kEdsCameraCommand_PressShutterButton, kEdsCameraCommand_ShutterButton_Completely);
			err = EdsSendCommand(camera, kEdsCameraCommand_PressShutterButton, kEdsCameraCommand_ShutterButton_OFF);
		}
		else if(c == 2){
			err = EdsSendCommand(camera, kEdsCameraCommand_PressShutterButton, kEdsCameraCommand_ShutterButton_Halfway);
			err = EdsSendCommand(camera, kEdsCameraCommand_PressShutterButton, kEdsCameraCommand_ShutterButton_OFF);
		}
		//else if(c == 3){
		//	err = EdsSendCommand(camera,
		//}
		else {
		
			
			
			//Sleep(300);
			//EdsSendCommand(camera,kEdsCameraCommand_TakePicture , 0);
			err = EdsSendCommand(camera, kEdsCameraCommand_PressShutterButton, kEdsCameraCommand_ShutterButton_Completely_NonAF);
			//Sleep(1200);
			err = EdsSendCommand(camera, kEdsCameraCommand_PressShutterButton, kEdsCameraCommand_ShutterButton_OFF);
			// Release camera
			
			
		
		}
	
	// Release camera
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
