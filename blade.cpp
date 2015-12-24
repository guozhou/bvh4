// bvh.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "Scene.h"
#include <tbb/task_scheduler_init.h>

int _tmain(int argc, _TCHAR* argv[])
{
	//tbb::task_scheduler_init init(tbb::task_scheduler_init::default_num_threads());
	try
	{
		Scene scene("CornellBox-Sphere.txt");
		scene.render();
		system((scene.fileName + ".exr").c_str());
	}
	catch (const std::string& err)
	{
		std::cerr << err << std::endl;
	}

	return 0;
}

