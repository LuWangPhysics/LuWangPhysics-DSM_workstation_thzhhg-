#ifndef PYTHON_PLOT_H
#define PYTHON_PLOT_H
//#include<Python.h>
struct python_plot{
   
    // PyObject* pModule;
    // PyObject*   pFunc;


    std::string save_path;
    
    
      python_plot(const std::string& save_path_ ){
          save_path=save_path_;
        //  Py_Initialize();
      //    PyRun_SimpleString("import sys");
         // PyRun_SimpleString("sys.path.append('/home/luwang/THz_DSM_code/DSM_THZ/my_plot')");
  

            
          // pModule = PyImport_ImportModule("read_plot");

          // //check whether the module is imported correctly
          //  if (pModule == NULL) {
          //     printf("ERROR importing module");
          //     PyErr_Print();
          //     exit(-1);
          //     }

          
    };
      ~ python_plot(){}
    // void plot_eff_spec(const std::string& b);
    // void plot_efield(const std::string& b);
    // void plot(const std::string& a,const std::string& b,const std::string& c);
     void clean();
};

// void python_plot::plot_eff_spec(const std::string& b){
//   pFunc = PyObject_GetAttrString(pModule, "plot1");

// //check whether the class is imported correctly
//  if (pFunc == NULL) {
//     printf("ERROR getting attribute");
//     exit(-1);
//     }
    
// //pass the saving string to python
//     PyObject_CallFunction(pFunc, "s,s",save_path.c_str(),b.c_str());
  
// };
    
// void python_plot::plot_efield(const std::string& b){
//    pFunc = PyObject_GetAttrString(pModule, "plot2");

//    // check whether the class is imported correctly
//      if (pFunc == NULL) {
//         printf("ERROR getting attribute");
//         exit(-1);
//         }
        
//    // pass the saving string to python
//    PyObject_CallFunction(pFunc, "s,s",save_path.c_str(), b.c_str());

// };


// void python_plot::plot(const std::string& a,const std::string& b,const std::string& c){
    
//   pFunc = PyObject_GetAttrString(pModule, "plot_general");

//      if (pFunc == NULL) {
//         printf("ERROR getting attribute");
//         exit(-1);
//         }
        

//     PyObject_CallFunction(pFunc, "s,s,s,s",save_path.c_str(),a.c_str(), b.c_str(),c.c_str());
 
// };
    
 void python_plot::clean(){
    
//        Py_DECREF(pFunc);
//        Py_DECREF(pModule);
//        Py_Finalize();
     };


    
    

#endif
