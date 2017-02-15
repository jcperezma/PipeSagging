#include <vector>
#include <list>

using namespace std;
template<class T> 
class row {
public:

	 virtual const T operator [](int i) const =0;    
    virtual T& operator[](int i)=0;
	virtual void resize(int n)=0;
	//virtual row<T>&  operator*(T& val)=0;		// times a scalar
	virtual void  rowPlus(row<T> & B, row<T> & result)=0;
	virtual void  rowMinus(row<T> & B, row<T> & result)=0;
	virtual T const  operator *(row<T> & vec)const=0;	// times a vector given as a row, very unlikely to happen
	virtual T const  operator*(vector<T>const & vec)const=0; //times a vector given as a vector
	virtual void  rowTimes(T & c, row<T> & result)=0; // times a scalar
	virtual T const  partialMultTo(vector<T>const & vec, int endColumn)const=0; //times a vector from column 0 to end Column 
	virtual T const  partialMultFrom(vector<T>const & vec, int initColumn)const=0; //times a vector from column initcolumn to the last 
	virtual void removeColumns(vector<int> & colIDs, vector<int> & idMap)=0;
	virtual void removeColumn(int colId)=0;
	virtual int const getSize()const =0;
	virtual void setToZero() =0;
protected:
};

template<class T> 
class Denserow: public row<T> {
	
public:
	Denserow(){};
	Denserow(int nRcolumns){rowData.resize(i)};
	~Denserow(){};

	/*
	Denserow(const Denserow<T> &B){
		rowData = B.getRowData();
	};*/

	Denserow<T> & operator = (Denserow<T> &B){
		rowData = B.getRowData();
		return *this;
	};

	void setToZero()override{
		for (int i = 0; i < rowData.size(); i++)
		{
			rowData[i] =0;
		}
};
	vector<T> const & getRowData(){
		return rowData;
	};
	const T operator[](int i) const override{
		return rowData[i];
	};
	T& operator[](int i)override{
		return rowData[i];
	};


	T const  operator*(vector<T>const & vec)const override{
		T result=0;
		int n = rowData.size();
		//cout<<rowData.size()<<endl;
		for (int i = 0; i < n; i++)
		{
			result += rowData[i]*vec[i];
			//cout <<result<<"  ";
		}
		return result;
	};
	
	void rowTimes(T & c,row<T> & result)override{

	result.resize(rowData.size());
	for (int i = 0; i < rowData.size(); i++)
		{
		result[i] = rowData[i] * c;
		}
	};


void rowPlus(row<T> & B,row<T> & result)override{
	// 

	result.resize(rowData.size());
	for (int i = 0; i < rowData.size(); i++)
	{
		result[i] = rowData[i] + B[i];
	}
	
}

void rowMinus(row<T> & B,row<T>&result)override{
	// returns a sparse row // constructing matrix C twice...

	result.resize(rowData.size());
	for (int i = 0; i < rowData.size(); i++)
	{
		result[i] = rowData[i] - B[i];
	}

}


    T const  partialMultTo(vector<T>const & vec, int  endColumn)const override{
		T result=0;
		for (int i = 0; i < endColumn; i++)
		{
			result += rowData[i]*vec[i];
		}
		return result;
	}; //times a vector from column 0 to end Column 
	T const  partialMultFrom(vector<T>const & vec, int  initColumn)const override {
		T result=0;
		for (int i = initColumn; i < rowData.size(); i++)
		{
			result += rowData[i]*vec[i];
		}
		return result;
	}; //times a vector from column initcolumn to the last 

	virtual T const  operator*(row<T> & vec)const override{
		T result=0;
		for (int i = 0; i < rowData.size(); i++)
		{
			result += rowData[i]*vec[i];
		}
		return result;
	};
	void resize(int n){
		rowData.resize(n);	
	};
	int const getSize()const{ return rowData.size(); };
	void removeColumns(vector<int> & colIDs, vector<int> & idMap) override { };
	void removeColumn(int colId){};
private:
	vector<T> rowData;
};

// Row element for Sparse Matrices
template<class T> 
class rowElement{
	T value;
	int column;
public:

	rowElement(const T& val=0, int col=-1): value(val),column(col){}; 
	rowElement(const rowElement&e): value(e.value),column(e.column){}; // copy constructor

	const rowElement& operator=(const rowElement&e){ 
		if(this != &e){
			value = e.value;
			column = e.column;
		}
		return *this;
	} // assignment operator


	~rowElement(){} // destructor
	
	T& getValue() { return value; };  // read the value
	const T getValue() const { return value; };  // read the value
	int getColumn() const{ return column; } // return the column
	const rowElement&  operator+=(const T&t){ // adding a T
		value += t;
		return *this;
	} 
	const rowElement& operator+=(const rowElement<T>&e){
		value += e.value;
		return *this;
	} // adding a rowElement
	const rowElement& operator-=(const T&t){
		value -= t;
		return *this;
	} // subtracting a T
	const rowElement& operator-=(const rowElement<T>&e){
		value -= e.value;
		return *this;
	} // subtracting a rowElement
	const rowElement& operator*=(const T&t){
		value *= t;
		return *this;
	} // multiplying by a T
	const rowElement& operator/=(const T&t){
		value /= t;
		return *this;
	} // dividing by a T
	const bool operator==(int col){
		return col ==column;
	}
	void setColumn(int newCol){
		column = newCol;
	}
	void setValue(T newVal){
		value = newVal;
	}
};


template<class T>
int operator<(const rowElement<T>&e, const rowElement<T>&f){
return e.getColumn() < f.getColumn();
} // smaller column index

template<class T>
int operator>(const rowElement<T>&e, const rowElement<T>&f){
return e.getColumn() > f.getColumn();
} // greater column index

template<class T>
int operator==(const rowElement<T>&e, const rowElement<T>&f){
return e.getColumn() == f.getColumn();
} // same column

template<class T>
const rowElement<T> operator+(const rowElement<T>&e, const T&t){
return rowElement<T>(e) += t;
} // rowElement plus a T

template<class T>
const rowElement<T> operator+(const T&t, const rowElement<T>&e){
return rowElement<T>(e) += t;
} // T plus rowElement

template<class T>
const rowElement<T> operator-(const rowElement<T>&e, const T&t){
return rowElement<T>(e) -= t;
} // rowElement minus T

template<class T>
const rowElement<T> operator*(const rowElement<T>&e, const T&t){
return rowElement<T>(e) *= t;
} // rowElement times a T
template<class T>
const rowElement<T> operator*(const T&t, const rowElement<T>&e){
return rowElement<T>(e) *= t;
} // T times rowElement
template<class T>
const rowElement<T> operator/(const rowElement<T>&e, const T&t){
return rowElement<T>(e) /= e;
} // rowElement divided by a T

template<typename T>
bool isZero (const rowElement<T>& element) { return (element.getValue()==0); }


template<typename T>
class IsColumn {
    int column;

public:
    IsColumn(int column) :
        column(column) {}

    bool operator()(const rowElement<T>& element) const {
		return column==element.getColumn();
    }
};


template<class T>
class Sparserow : public row<T>{
public:

Sparserow(const T&val=-1,int col=-1){
	rowData.push_back(rowElement<T>(val,col));
	
} // constructor

/*
Sparserow( Sparserow<T> & spRow){
	spRow.getRowData(rowData);
	
}
*/

/*
const Sparserow<T>& operator=( Sparserow<T>&e){ 
		if(this != &e){	
	for (int i = 0; i < e.getSize(); i++)
	{
		insertNextItem(e[i],3 );
	}
		}
		return *this;
	} // assignment operator
	*/
void setToZero()override{
	for ( auto & e:rowData)
	{
		e.setValue(0);
	}
};
T&  operator[](int i) override{
for ( auto it =rowData.begin(); it !=rowData.end(); ++it)
	{
		
		int column = it->getColumn();
		if(it->getColumn()==i) 	
			return it->getValue();
			//return it->getValue();
	}
// didn't find column

	T result = 0; // this is terribly wrong
	return result;
}

void rowPlus(row<T> & B, row<T> &result) override{
	// returns a sparse row
	
	
	for (auto & e:rowData)
	{
		if (e.getColumn()!=-1)
			result[e.getColumn()] = e.getValue() + B[e.getColumn()];
	}

}




void  rowTimes(T & c, row<T> & result)override{
	// returns a sparse row
	// result has been initialized prior to this call
	
	for (auto & e:rowData)
	{
		if (e.getColumn()!=-1)
			result[e.getColumn()] = e.getValue() *c;
	}

}

void rowMinus(row<T> & B,row<T> & result)override{
	// returns a sparse row
	
	for (auto & e:rowData)
	{
		if (e.getColumn()!=-1)
			result[e.getColumn()] = e.getValue() - B[e.getColumn()];
	}


}


const T  operator[](int i) const override{
for ( auto it =rowData.begin(); it !=rowData.end(); ++it)
	{
		
		int column = it->getColumn();
		if(it->getColumn()==i) 	
			return it->getValue();
			//return it->getValue();
	}
// didn't find column

	T result =0;
	return result;
}

rowElement<T> & getElementAt(int i ) {
	auto it= rowData.begin();
	advance(it,i);
	 return *it;
}

 void removeZeros(){
	 rowData.remove_if(isZero<T>);
 };

 void removeColumns(vector<int> & colIDs, vector<int> & idMap) override{
	 for (auto colID:colIDs )
	 {
		 
		 rowData.remove_if(IsColumn<T>(colID));
	 }

	 // update columns Ids

	 for (auto & e: rowData)
	 {
		 if (e.getColumn()!=-1)
			e.setColumn(e.getColumn()- idMap[e.getColumn()]  );
	 }

 }

 void removeColumn(int column) override{
	 rowData.remove_if(IsColumn<T>(column));
 }

 void getRowData(list<rowElement<T>> & rowD) {
	rowD.insert(rowD.end(), rowData.begin(), rowData.end());
}
void resize(int n){
		rowData.resize(n);	
	};


void insertNextItem(const T&val, int col){
	rowElement<T> e(val,col);
	if (e<*rowData.begin())
	{
		rowData.emplace_front(e);
		return;
	}
	for ( auto it =rowData.begin(); it !=rowData.end(); ++it)
	{
		if (*it==e) // there is an element in that position, give it a new value
		{
			*it = rowElement<T>(val,col);
		}
		if (it!=prev(rowData.end()))
		{
			if (*it < e && e < *next(it))
			{
				rowData.emplace(++it,e);
				
				return;
			}
		} else
		{
			if (*it < e )
			{
				rowData.emplace(++it,e);
				
				return;
			}
		}
		
	}
} 

T const  operator*(vector<T>const & vec)const override{
	T result =0;
	for ( auto e: rowData)
	{
		if (e.getColumn()>=0)
		{
			result+=e.getValue()*vec[e.getColumn()];	
		}	
	}
	return result;
}
T const  partialMultTo(vector<T>const & vec, int  endColumn)const override{
	T result =0;
	for ( auto e: rowData)
	{
		if (e.getColumn()>=0 && e.getColumn()<endColumn )
		{
			result+=e.getValue()*vec[e.getColumn()];	
		}	
	}
	return result;
}; //times a vector from column 0 to end Column 

T const  partialMultFrom(vector<T>const & vec, int  initColumn)const override{
	T result =0;
	for ( auto e: rowData)
	{
		if (e.getColumn()>=initColumn  )
		{
			result+=e.getValue()*vec[e.getColumn()];	
		}	
	}
	return result;
}; //times a vector from column initcolumn to the last 
	

T const   operator*(row<T> & vec)const override{
		T result =0;
	for ( auto e: rowData)
	{
		result+=e.getValue()*vec[e.getColumn()];
	}
	return result;
};



const T rowSum() const{
	T result =0;
	for (auto e : rowData)
	{
		result+=e.getValue();
	}
return result;
	
}

int const getSize()const{ return rowData.size(); };



protected:
	list<rowElement<T>> rowData;
	
};

/* revise these
template<class T>
const row<T> operator*(const row<T>&r, const T&t){
return row<T>(r) *= t;
} // row times T

template<class T>
const row<T> operator*(const T&t, const row<T>&r){
return row<T>(r) *= t;
} // T times row

template<class T>
const row<T> operator/(const row<T>&r, const T&t){
return row<T>(r) /= t;
} // row divided by a T

*/
template<typename T>
Sparserow<T>  operator +(Sparserow<T> & A , Sparserow<T> & B) {
	// returns a sparse row
	Sparserow<T>C = *this;
	
	for (auto & e:rowData)
	{
		if (e.getColumn()!=-1)
			C[e.getColumn()] = e.getValue() + B[e.getColumn()];
	}

	return C;
}