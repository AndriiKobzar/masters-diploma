using System;
using System.Net.Http;
using System.Threading.Tasks;
using System.Web;

namespace DiplomaProject
{
    public class WolframClient
    {
        const string API_KEY = "KVQX75-XGT5TGH3R8";
        
        
        public async Task ExecuteRequest(string request)
        {
            using (var httpClient = new HttpClient())
            {
                using (var httpRequest = await httpClient.GetAsync(GetRequestString(request)))
                {
                    httpRequest.EnsureSuccessStatusCode();
                    Console.WriteLine(await httpRequest.Content.ReadAsStringAsync());
                }
            }
        }

        private string GetRequestString(string request)
        {
            string encodedRequest = HttpUtility.UrlEncode(request);
            return $"http://api.wolframalpha.com/v2/query?input={encodedRequest}&appid={API_KEY}";
        }
    }
}