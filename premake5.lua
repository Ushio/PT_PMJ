include "libs/PrLib"

workspace "SimpleRT"
    location "build"
    configurations { "Debug", "Release" }
    startproject "main"

architecture "x86_64"

externalproject "prlib"
	location "libs/PrLib/build" 
    kind "StaticLib"
    language "C++"

project "main"
    kind "ConsoleApp"
    language "C++"
    targetdir "bin/"
    systemversion "latest"
    flags { "MultiProcessorCompile", "NoPCH" }

    -- Src
    files { "main.cpp", "EzEmbree.hpp" }

    -- rapidjson
    includedirs { "libs/rapidjson/include" }
    files { "libs/rapidjson/include/**.h" }

    -- embree
    includedirs { "libs/embree/include" }
    libdirs { "libs/embree/lib" }
    filter {"Debug"}
        links { "embree3", "tbb_debug", "tbbmalloc_debug" }
    filter {"Release"}
        links { "embree3",  "tbb", "tbbmalloc" }
    filter{}
    postbuildcommands { 
        "{COPY} %{prj.location}../libs/embree/bin/embree3.dll %{cfg.targetdir}",
        "{COPY} %{prj.location}../libs/embree/bin/tbb.dll %{cfg.targetdir}",
        "{COPY} %{prj.location}../libs/embree/bin/tbb_debug.dll %{cfg.targetdir}",
        "{COPY} %{prj.location}../libs/embree/bin/tbbmalloc.dll %{cfg.targetdir}",
        "{COPY} %{prj.location}../libs/embree/bin/tbbmalloc_debug.dll %{cfg.targetdir}",
     }

    -- prlib
    -- setup command
    -- git submodule add https://github.com/Ushio/prlib libs/prlib
    -- premake5 vs2017
    dependson { "prlib" }
    includedirs { "libs/prlib/src" }
    libdirs { "libs/prlib/bin" }
    filter {"Debug"}
        links { "prlib_d" }
    filter {"Release"}
        links { "prlib" }
    filter{}

    symbols "On"

    filter {"Debug"}
        runtime "Debug"
        targetname ("Main_Debug")
        optimize "Off"
    filter {"Release"}
        runtime "Release"
        targetname ("Main")
        optimize "Full"
    filter{}
