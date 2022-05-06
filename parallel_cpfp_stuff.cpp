#include "parallel_cpfp_stuff.h"

#define WRITE_IN
  
#ifdef USING_MPI
	int EASY_SEND_AND_RECEIVE::receive_blocking_signal_as_tag_and_number(int source){
//         IFMASTER cout << "ss " << __LINE__ << endl;
  
     	v_request.resize(2,MPI_REQUEST_NULL);
//         IFMASTER cout << "ss " << __LINE__ << endl;
        
		v_status.resize(2);
//         IFMASTER cout << "ss " << __LINE__ << endl;
        MPI_Recv(&send_recv_tag,1,MPI_INT,source,TAG_SENDING_SIGNAL, MPI_COMM_WORLD, &v_status[0]);
//         IFMASTER cout << "ss " << __LINE__ << endl;
        MPI_Recv(&number_tag,1,MPI_INT,source,TAG_SENDING_SIGNAL, MPI_COMM_WORLD, &v_status[1]);

        //EASY_Recv(&send_recv_tag,1,EASY_INT,source,TAG_SENDING_SIGNAL);
//         IFMASTER cout << "ss " << __LINE__ << endl;
        //EASY_Recv(&number_tag,1,EASY_INT,source,TAG_SENDING_SIGNAL);

        double_tag = 0;
//         IFMASTER cout << "ss " << __LINE__ << endl;
        
		return 0;
        
	}

	int EASY_SEND_AND_RECEIVE::receive_blocking_signal_as_tag_and_double(int source){
		v_request.resize(2,MPI_REQUEST_NULL);
		v_status.resize(2);
		MPI_Recv(&send_recv_tag,1,MPI_INT,source,TAG_SENDING_SIGNAL, MPI_COMM_WORLD, &v_status[0]);
		MPI_Recv(&double_tag,1,MPI_DOUBLE,source,TAG_SENDING_SIGNAL, MPI_COMM_WORLD, &v_status[1]);
		number_tag = 0;
		return 0;
	}

	int EASY_SEND_AND_RECEIVE::send_blocking_signal_as_tag_and_number(int tag, int number, int destination){
		send_recv_tag = tag;
		number_tag = number;
		v_request.resize(2,MPI_REQUEST_NULL);
		v_status.resize(2);
		MPI_Send(&send_recv_tag,1,MPI_INT,destination,TAG_SENDING_SIGNAL, MPI_COMM_WORLD);
		MPI_Send(&number_tag,1,MPI_INT,destination,TAG_SENDING_SIGNAL, MPI_COMM_WORLD);
		double_tag = 0;
		return 0;
	}

	int EASY_SEND_AND_RECEIVE::send_blocking_signal_as_tag_and_double(int tag, double d_number, int destination){
		send_recv_tag = tag;
		double_tag = d_number;
		v_request.resize(2,MPI_REQUEST_NULL);
		v_status.resize(2);
		MPI_Send(&send_recv_tag,1,MPI_INT,destination,TAG_SENDING_SIGNAL, MPI_COMM_WORLD);
		MPI_Send(&double_tag,1,MPI_DOUBLE,destination,TAG_SENDING_SIGNAL, MPI_COMM_WORLD);
		number_tag = 0;
		return 0;

	}


	int EASY_SEND_AND_RECEIVE::get_tag() const{ return send_recv_tag;}
	int EASY_SEND_AND_RECEIVE::get_number() const {return number_tag;}
	double EASY_SEND_AND_RECEIVE::get_double() const {return double_tag;}

void EASY_SEND_AND_RECEIVE::print(ostream& out) const{
	out << "Tag: " << send_recv_tag ;
	out << " number_tag: " << number_tag;
	out << " double_tag: " << double_tag << endl;
}


#else
 double my_signaler::get_double() const {return double_tag;}
 int my_signaler::get_number() const { return number_tag;}
 int my_signaler::get_tag() const {return send_recv_tag;}
 int my_signaler::get_MY_WORLD_SIZE() const {return MY_WORLD_SIZE;}
 int my_signaler::get_MY_SLAVE_ID() const {return MY_SLAVE_ID;}
 
 void my_signaler::set_tag(int tag)
{
 send_recv_tag = tag;
 return;
}

 
 void my_signaler::print(std::ostream& out) const{
    out << "MY_SLAVE_ID:" << MY_SLAVE_ID;
	out << "Tag: " << send_recv_tag ;
	out << " number_tag: " << number_tag;
	out << " double_tag: " << double_tag << endl;
 
}



int my_signaler::object_check(double _ms)
{
 if(to.pause()>_ms){
  to.reset();
 }
 else{
  return -1;
 }
 
    string fname = object_receiver_filename.str();
    ifstream myfile (object_receiver_filename.str(),ios::in);
    if (!file_exists(fname)) return -1;
    string dummy; 
    //cout << " printing at " << __FILE__ << " and line "<< __LINE__ <<endl; 
    int number_of_new_objects  = 0; 
    for (int i = 0; i < object_list.size(); i++){
     //cout << " printing at " << __FILE__ << " and line "<< __LINE__ <<endl; 
         std::getline (myfile,dummy);
     
    }
    //cout << " printing at " << __FILE__ << " and line "<< __LINE__ <<endl; 
    string tsource, ttag, tval;
    //cout << " printing at " << __FILE__ << " and line "<< __LINE__ <<endl; 
    while (!myfile.eof()){
     //cout << " printing at " << __FILE__ << " and line "<< __LINE__ <<endl; 
       std::getline (myfile,dummy);
       // << " printing at " << __FILE__ << " and line "<< __LINE__ << " dummy " << dummy <<endl; 
       stringstream s; 
       //cout << " printing at " << __FILE__ << " and line "<< __LINE__ <<endl; 
       
       s<< dummy;
       //cout << " printing at " << __FILE__ << " and line "<< __LINE__ <<endl; 
       s  >> tsource >> ttag;
       //cout << " printing at " << __FILE__ << " and line "<< __LINE__ << endl; 
       object ob (stoi(tsource),stoi(ttag),s.str());
       //cout << " printing at " << __FILE__ << " and line "<< __LINE__ <<endl; 
       object_list.push_back(ob);
       number_of_new_objects++;
    }
      
    myfile.close();
    return number_of_new_objects;
}


int my_signaler::signal_check(double _ms){ 
//     static int isss = 0;
//     isss++;
    if(ts.pause()>_ms){
     
       ts.reset();
    }
    else{
     //if (isss>10) return TAG_SIGNALER_NO_NEW_SIGNAL_TIME_LIMIT;
//      cout << ts.pause() << " " << _ms << endl;
//      cout << "check_signal returns TAG_SIGNALER_NO_NEW_SIGNAL_TIME_LIMIT" << TAG_SIGNALER_NO_NEW_SIGNAL_TIME_LIMIT << endl;
     return TAG_SIGNALER_NO_NEW_SIGNAL_TIME_LIMIT;
    }
    
    string filename = signal_receiver_filename.str();
    if (!file_exists(filename)){
      //cout << "file "  << signal_receiver_filename.str() << " does not exist my_signaler::signal_check returns -1" <<endl; 
     return -1;
    }
    ifstream myfile (signal_receiver_filename.str(),ios::in);
    string dummy; 
    int number_of_new_signals  = 0; 
    for (int i = 0; i < signal_list.size(); i++){
         
         std::getline (myfile,dummy);
         //cout << __LINE__ << __FILE__ << "i: " <<i << endl;
    }
    
    string tsource, ttag, tval;
    while (!myfile.eof()){
       std::getline (myfile,dummy);
       if (dummy == "") break;
       stringstream s; 
       s<< dummy;
//       cout << __LINE__ << __FILE__ << "n: " <<number_of_new_signals << endl;
       s >> tsource >> ttag >> tval;
//        cout << __LINE__ << __FILE__ << "s: " <<MY_SLAVE_ID << endl;
       
//        cout << __LINE__ << __FILE__ << "d: " <<dummy << endl;
       
//        cout << __LINE__ << __FILE__ <<  ttag << tval << tsource << endl;
       
       if (stoi(ttag)>DOUBLE_TAG_LIMIT){
           double dd = -1;
           stringstream ss(tval);
           ss >> dd;
           signal sd(stoi(tsource),stoi(ttag),stod(tval));
           signal_list.push_back(sd);
//            cout << "printing sd" <<endl;
//            sd.print();
//                   cout << __LINE__ << __FILE__ << "s: " <<MY_SLAVE_ID << endl;
//            
//                   cout << __LINE__ << __FILE__ << "d: " <<dummy << endl;
//            
//                   cout << __LINE__ << __FILE__ <<  ttag << " tval:"<< tval<< " " << tsource << " dd " << dd <<endl;
//            
//            cout << "printing sd" <<endl;
           
        
       }
       else {
           signal si(stoi(tsource),stoi(ttag),stoi(tval));
           //si.print();
           signal_list.push_back(si);
        
       }
       
//        cout << __LINE__ << __FILE__ << endl;
       number_of_new_signals++;
    }
//     cout << __LINE__ << __FILE__ << "n: " <<number_of_new_signals << endl;
    myfile.close();
//     cout << "check_signal returns " << number_of_new_signals << endl;
    return number_of_new_signals;
    
};

int my_signaler::receive_signal_as_tag_and_number(int source, int &flag, ALT1_RECV_TYPE RECV_TYPE )
{   
    flag = TAG_SIGNALER_RETURNS_NOERROR;
    int iso = signal_check();
 
    if (source == -1){
        for (int i = int(signal_list.size()) - 1; i>=0; --i){
           if (signal_list[i].intmi){
             signal_list[i].used = true;
             return signal_list[i].intval;
           }
        }
    }
    else{
     for (int i = int(signal_list.size()) - 1; i>=0; --i){
           if (signal_list[i].intmi && signal_list[i].source == source){
             signal_list[i].used = true;
             return signal_list[i].intval;
           }
        }
    }
    
    flag = TAG_SIGNALER_RETURNS_ERROR;
    return TAG_SIGNALER_RETURNS_ERROR;
    
 
}


double my_signaler::receive_signal_as_tag_and_double(int source, int &flag, ALT1_RECV_TYPE RECV_TYPE)
{
    flag = TAG_SIGNALER_RETURNS_NOERROR;
    
    int iso = signal_check();
    if (source == -1){
        for (int i = int(signal_list.size()) - 1; i>=0; --i){
           if (!signal_list[i].intmi){
             signal_list[i].used = true;
             return signal_list[i].doubleval;
           }
        }
    }
    else{
     for (int i = int(signal_list.size()) - 1; i>=0; --i){
           if (!signal_list[i].intmi && signal_list[i].source == source){
             signal_list[i].used = true;
             return signal_list[i].doubleval;
           }
        }
    }
    flag = TAG_SIGNALER_RETURNS_ERROR;
    return TAG_SIGNALER_RETURNS_ERROR;
    
 
}

std::string my_signaler::recv_object(int source, int tag, ALT1_RECV_TYPE RECV_TYPE)
{
    int iso = object_check();
    if (source == -1){
        for (int i = int(object_list.size()) - 1; i>=0; --i){
         if (object_list[i].tag == tag && ! (object_list[i].used) ){
             object_list[i].used = true;
             return object_list[i].object_string;
           }
        }
    }
    
    else{
     for (int i = int(signal_list.size()) - 1; i>=0; --i){
      if (object_list[i].tag == tag && ! (object_list[i].used) && object_list[i].source == source){
             object_list[i].used = true;
             return object_list[i].object_string;
           }
        }
    }
    return "";
}


my_signaler::my_signaler()
{
 
}

my_signaler::my_signaler(const my_signaler& other)
{
 ts.start();
 to.start();
 communication_world_uniq_file_prefix = other.communication_world_uniq_file_prefix;
 MY_SLAVE_ID = other.MY_SLAVE_ID;
 MY_WORLD_SIZE = other.MY_WORLD_SIZE;
  signal_receiver_filename << communication_world_uniq_file_prefix << "_SLAVE" << to_string(MY_SLAVE_ID) << "_signal.in" ;
 signal_sender_filename << communication_world_uniq_file_prefix << "_SLAVE" << to_string(MY_SLAVE_ID) << "_signal.out" ;
 object_receiver_filename << communication_world_uniq_file_prefix << "_SLAVE" << to_string(MY_SLAVE_ID) << "_object.in" ;
 object_sender_filename << communication_world_uniq_file_prefix << "_SLAVE" << to_string(MY_SLAVE_ID) << "_object.out" ;

}


my_signaler &my_signaler::operator= (const my_signaler& other)
{
 ts.start();
 to.start();
 communication_world_uniq_file_prefix = other.communication_world_uniq_file_prefix;
 MY_SLAVE_ID = other.MY_SLAVE_ID;
 MY_WORLD_SIZE = other.MY_WORLD_SIZE;
  signal_receiver_filename << communication_world_uniq_file_prefix << "_SLAVE" << to_string(MY_SLAVE_ID) << "_signal.in" ;
 signal_sender_filename << communication_world_uniq_file_prefix << "_SLAVE" << to_string(MY_SLAVE_ID) << "_signal.out" ;
 object_receiver_filename << communication_world_uniq_file_prefix << "_SLAVE" << to_string(MY_SLAVE_ID) << "_object.in" ;
 object_sender_filename << communication_world_uniq_file_prefix << "_SLAVE" << to_string(MY_SLAVE_ID) << "_object.out" ;
return *this;
}
my_signaler::my_signaler(std::string s, int ID, int WS)
{
 ts.start();
 to.start();
 communication_world_uniq_file_prefix = s;
 MY_SLAVE_ID = ID;
 MY_WORLD_SIZE = WS;
 signal_receiver_filename << communication_world_uniq_file_prefix << "_SLAVE" << to_string(MY_SLAVE_ID) << "_signal.in" ;
 signal_sender_filename << communication_world_uniq_file_prefix << "_SLAVE" << to_string(MY_SLAVE_ID) << "_signal.out" ;
 object_receiver_filename << communication_world_uniq_file_prefix << "_SLAVE" << to_string(MY_SLAVE_ID) << "_object.in" ;
 object_sender_filename << communication_world_uniq_file_prefix << "_SLAVE" << to_string(MY_SLAVE_ID) << "_object.out" ;
 
 
}


string my_signaler::generate_slaveX_signal_in(int sl_id){
 stringstream s;
 s << communication_world_uniq_file_prefix << "_SLAVE" << to_string(sl_id) << "_signal.in" ;
 return s.str();
}

string my_signaler::generate_slaveX_object_in(int sl_id)
{
 stringstream s;
 s << communication_world_uniq_file_prefix << "_SLAVE" << to_string(sl_id) << "_object.in" ;
 return s.str();
}


int my_signaler::receive_signal_i(int source, int tag, int &flag, ALT1_RECV_TYPE RECV_TYPE) // returns -1 if for used messages 
{
    flag = TAG_SIGNALER_RETURNS_NOERROR;
     int iso = signal_check();
    if (source == -1){
        for (int i = int(signal_list.size()) - 1; i>=0; --i){
           if (signal_list[i].intmi && !signal_list[i].used && signal_list[i].tag == tag){
             signal_list[i].used = true;
             return signal_list[i].intval;
           }
        }
    }
    else{
     for (int i = int(signal_list.size()) - 1; i>=0; --i){
           if (signal_list[i].intmi && !signal_list[i].used && signal_list[i].tag == tag && signal_list[i].source == source){
             signal_list[i].used = true;
             return signal_list[i].intval;
           }
        }
    }
    flag = TAG_SIGNALER_RETURNS_ERROR;
    return flag;
 
}

void signal::print(){
 cout << "source: " << source << " tag: "<< tag << " intval: " <<intval  << " doubleval: " << doubleval << " used: " <<used<< " intmi: "<< intmi <<endl;   
}
double my_signaler::receive_signal_d(int source, int tag, int &flag, ALT1_RECV_TYPE RECV_TYPE) // returns -1 if for used messages 
{
    //blocking calismiyor/
 
    flag = TAG_SIGNALER_NEW_SIGNAL;
    int iso = 0;
    if (RECV_TYPE == ALT1_RECV_BLOCKING ){
        while (1){
         if (source == -1){
          for (int i = int(signal_list.size()) - 1; i>=0; --i){
           //cout << "signal i:"<< i << " with tag " <<   signal_list[i].tag << " and source " <<  signal_list[i].source  << " is used: " <<  signal_list[i].used << endl;
           if (!(signal_list[i].intmi) && !(signal_list[i].used) && signal_list[i].tag == tag){
            signal_list[i].used = true;
            return signal_list[i].doubleval;
           }
          }
          
         }
         else{
          for (int i = int(signal_list.size()) - 1; i>=0; --i){
           
           //cout << "signal i:"<< i << " with tag " <<   signal_list[i].tag << " and source " <<  signal_list[i].source  << " is used: " <<  signal_list[i].used << endl;
           //signal_list[i].print();
           if (!(signal_list[i].intmi) && !(signal_list[i].used) && signal_list[i].tag == tag && signal_list[i].source == source){
            signal_list[i].used = true;
            
            return signal_list[i].doubleval;
           }
          }
         }
         
         while (iso<=0) {
             iso = signal_check(0);

             //cout << "iso" << iso << endl;
         }
         //cout << "iso" << iso << endl;
        }
    }

    iso = signal_check();
    
    //cout << __LINE__ << " iso " << iso << endl;
    if (source == -1){
        for (int i = int(signal_list.size()) - 1; i>=0; --i){
         //cout << "signal i:"<< i << " with tag " <<   signal_list[i].tag << " and source " <<  signal_list[i].source  << " is used: " <<  signal_list[i].used << endl;
           if (!(signal_list[i].intmi) && !(signal_list[i].used) && signal_list[i].tag == tag){
             signal_list[i].used = true;
             return signal_list[i].doubleval;
           }
        }
        
    }
    else{
     for (int i = int(signal_list.size()) - 1; i>=0; --i){
      //cout << "signal i:"<< i << " with tag " <<   signal_list[i].tag << " and source " <<  signal_list[i].source  << " is used: " <<  signal_list[i].used << endl;
      
      //signal_list[i].print();
      if (!(signal_list[i].intmi) && !(signal_list[i].used) && signal_list[i].tag == tag && signal_list[i].source == source){
             signal_list[i].used = true;
             
             return signal_list[i].doubleval;
           }
        }
    }
    
    flag = TAG_SIGNALER_RETURNS_ERROR;
    return flag;
 
}

int my_signaler::send_object(int destination, int tag, std::string objs, ALT1_SEND_TYPE SEND_TYPE)
{
 stringstream o;
 o << destination << " " << tag << " " << objs;
 ofstream myfile (object_sender_filename.str(),ios::app);
 myfile << o.str() << endl;
 myfile.close();
 /* create a copy of input file of the destination here 
  */
 #ifdef WRITE_IN
 stringstream s;
 s << MY_SLAVE_ID << " " << tag << " " << objs;

 if (destination == DESTINATION_ALL_SLAVES){
  for (int i =1; i < MY_WORLD_SIZE; ++i){
   string s2 = generate_slaveX_object_in(i);
   ofstream myfile2 (s2,ios::app);
   myfile2 << s.str() << endl;
   myfile2.close();
  }
 }
 else{
  string s3 = generate_slaveX_object_in(destination);
  ofstream myfile3 (s3,ios::app);
  myfile3 << s.str() << endl;
  myfile3.close();
 }
#endif
 
 return 1;
 
}

int my_signaler::send_signal_as_tag_and_double(int destination, int tag, double d_number)
{
  stringstream o;
  o << destination << " " << tag << " " << d_number;
 ofstream myfile (signal_sender_filename.str(),ios::app);
 myfile << o.str() << endl;
 myfile.close();
 
 /* create a copy of input file of the destination here 
  */
 #ifdef WRITE_IN 
 stringstream s;
 s << MY_SLAVE_ID << " " << tag << " " << d_number;
 
 if (destination == DESTINATION_ALL_SLAVES){
  for (int i =1; i < MY_WORLD_SIZE; ++i){
   string s2 = generate_slaveX_signal_in(i);
   ofstream myfile2 (s2,ios::app);
   myfile2 << s.str() << endl;
   myfile2.close();
  }
 }
 else{
  string s3 = generate_slaveX_signal_in(destination);
  ofstream myfile3 (s3,ios::app);
  myfile3 << s.str() << endl;
  myfile3.close();
 }
#endif
  return 1;
}

void my_signaler::clean_files(bool sig, bool obj)
{
 
 if (sig) {
  unlink(signal_sender_filename.str().c_str());
#ifdef WRITE_IN
  unlink(signal_receiver_filename.str().c_str());
#endif
  
 }
 if (obj) {
  unlink(object_sender_filename.str().c_str());

#ifdef WRITE_IN
  unlink(object_receiver_filename.str().c_str());
#endif

 }
 return;

}


int my_signaler::send_signal_as_tag_and_number(int destination, int tag, int number)
{
  stringstream o;
  o << destination << " " << tag << " " << number;
 ofstream myfile (signal_sender_filename.str(),ios::app);
 myfile << o.str() << endl;
 myfile.close();
 
 
 /* create a copy of input file of the destination here 
  */
 #ifdef WRITE_IN
 stringstream s;
 s << MY_SLAVE_ID << " " << tag << " " << number;
 
 if (destination == DESTINATION_ALL_SLAVES){
  for (int i =1; i < MY_WORLD_SIZE; ++i){
   string s2 = generate_slaveX_signal_in(i);
   ofstream myfile2 (s2,ios::app);
   myfile2 << s.str() << endl;
   myfile2.close();
  }
 }
 else{
  string s3 = generate_slaveX_signal_in(destination);
  ofstream myfile3 (s3,ios::app);
  myfile3 << s.str() << endl;
  myfile3.close();
 }
#endif
  return 1;
}

int my_signaler::check_signal(int source, int tag)
{
 int ois = signal_check();
 
 if (source == -1){
  for (int i = int(signal_list.size()) - 1; i>=0; --i){
   if (!(signal_list[i].used) && signal_list[i].tag == tag){
    //std::cout << "signal i:"<< i << " with tag " <<   signal_list[i].tag << " and source " <<  signal_list[i].source  << " is used: " <<  signal_list[i].used << endl;
    return TAG_SIGNALER_PROBE_SUCCESSFUL;
   }
  }
  
 }
 else{
  for (int i = int(signal_list.size()) - 1; i>=0; --i){
   if (!(signal_list[i].used) && signal_list[i].tag == tag && signal_list[i].source == source){
    //std::cout << "signal i:"<< i << " with tag " <<   signal_list[i].tag << " and source " <<  signal_list[i].source  << " is used: " <<  signal_list[i].used << endl;
    return TAG_SIGNALER_PROBE_SUCCESSFUL;
   }
  }
 }

 
 return TAG_SIGNALER_NO_NEW_SIGNAL;
 
}

/*
		int  send_signal_as_tag_and_number(int destination, int tag, int number); 
		int  send_signal_as_tag_and_double(int destination, int tag, double d_number);
        
	    int check_signal(int source, int tag); //returns 1 if there is a new signal tag -1 for any signal 
        int receive_signal_i(int source, int tag);
        double receive_signal_d(int source, int tag);
}

*/


object::object(){}
object::object(int s, int t, std::string ss)
{
  source = s;
  tag = t;
  object_string = ss;
  used = false;
 
}
signal::signal()
{
}

signal::signal(int s, int t, double d)
{
 source = s;
 tag = t;
 intmi = false;  
 doubleval = d;
 intval = -999;
 used = false;
 
}

signal::signal(int s, int t, int i)
{
  source = s;
 tag = t;
 intmi = true;  
 doubleval = -999;
 intval = i;
 used = false;
 
}


#endif


// kate: indent-mode cstyle; indent-width 1; replace-tabs on; 
