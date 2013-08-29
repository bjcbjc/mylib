function obj = DAVIDWebService

obj.endpoint = 'http://david.abcc.ncifcrf.gov:80/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap11Endpoint/';
obj.wsdl = 'http://david.abcc.ncifcrf.gov/webservice/services/DAVIDWebService?wsdl';

obj = class(obj,'DAVIDWebService');

