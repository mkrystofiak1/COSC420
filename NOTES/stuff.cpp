#include<iostream>

#ifdef __cplusplus
class A {
	private:
		int num;
		std::string name[255];
		double price;
	
	public:
		A(int n, std::string name, double price): num(n), name(name), price(price) {};

		std::string toString(){
			return name + " " + std::to_string(num) + " " + std::to_string(price);
		}

};

#else
typdef struct A;

void makeA(A* obj){
}
extern "C" void f() {
	std::cout << "Hello from function f!" << std::endl;
}


